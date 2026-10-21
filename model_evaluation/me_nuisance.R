##########
# title: nuisance-evaluation pipelines - model evaluation study
##########
# A second, independent nuisance-estimation pipeline used only to score the 9
# candidate models from me_models.R against each other and against true PEHE
# - not to fit them. No other study in this repo needs this concept, so it
# doesn't map onto config/dgms/models/analysis/check/collect/metrics and gets
# its own file, the same way crossfitting/ has files beyond that 7-role
# floor (cf_testing.R, cf_results.R, ...)
#
# Two independent estimators of the same 2 nuisance targets (mu_DR, pi):
# XGBoost with a hand-tuned CV grid, and H2O AutoML, each under the arms
# defined in me_config.R's NUISANCE_ARMS. cv_shared reuses the candidates' own
# folds (ADEMP.md, "Why there is no cv_indep arm").
#
# `holdout` fits and predicts on the fold the candidate held out: not row-level
# honest, but fully decoupled from the candidate. One small block cannot give
# both properties.
#
# Outcome models are per arm (T-learner): mu0_DR on the W = 0 rows and mu1_DR on
# the W = 1 rows of each training set, in both pipelines and every arm, matching
# the candidates (DEVLOG bug R).
require(h2o)
require(xgboost)

unlist_order <- function(reslist) {
  df <- do.call(rbind, reslist)
  df <- df %>%
    arrange(row_index) %>%
    select(-row_index)
  return(df)
}

#' The AIPW pseudo-outcomes from one arm's nuisance predictions
#'
#' Propensity trimmed to [0.05, 0.95] (trim_ps(), R/utils.R), as the 9
#' candidates' and every other study's are. Untrimmed, AutoML's propensity can
#' hit exactly 0 or 1, making phi NaN or +-Inf (DEVLOG bug S). The truth is
#' 0.5, so trimming only moves an estimate that is already wrong toward it. `pi` keeps the raw prediction (me_strategies.R's
#' pi_diagnostics read it); `pi_trim` is what phi and the IF score use.
#'
#' Row-wise and idempotent on its own output: `pi` is never overwritten, so
#' me_metrics.R recomputes phi from the saved mu0_DR/mu1_DR/pi with this same
#' function, and results written before the trimming score as if fit after.
calculate_pseudos <- function(df, Y, W) {
  df <- df %>%
    mutate(
      pi_trim = trim_ps(pi),
      mu_DR = mu1_DR * W + (1 - W) * mu0_DR
    )
  df <- df %>%
    mutate(
      phi = mu1_DR - mu0_DR +
        ((Y - mu_DR) * (W - pi_trim)) / (pi_trim * (1 - pi_trim)),
      # phi05 fixes the propensity at 0.5 instead of using the estimated pi.
      # It is not an approximation: R/dgm_scenarios.R assigns treatment with
      # `W <- rbinom(n, 1, 0.5)`, so 0.5 IS the true propensity and phi05 is
      # the oracle-pi pseudo-outcome. me_metrics.R's dr05_* and cal05q*_*
      # score columns are built from it, and the dr/dr05 contrast is what
      # isolates the cost of estimating a propensity that did not need
      # estimating. Computed here since this study's first commit - which is
      # why those columns needed no rerun to add.
      phi05 = mu1_DR - mu0_DR + ((Y - mu_DR) * (W - 0.5)) / (0.5 * (1 - 0.5))
    )
  return(df)
}

# XGB
create_param_grid <- function(n_cores = 1) {
  grid_base <- expand.grid(
    nthread = n_cores,
    eta = c(0.1, 0.2, 0.3),
    max_depth = c(2, 4, 6),
    subsample = c(0.6, 0.8),
    colsample_bytree = c(0.7, 1)
  )

  list(
    continuous = grid_base %>% mutate(objective = "reg:squarederror"),
    binary = grid_base %>% mutate(objective = "binary:logistic")
  )
}

# Fewest training rows either treatment arm may have before a per-arm outcome
# fit is refused. XGBoost's inner CV is 5-fold, so below this each inner
# fold trains on a handful of rows. The thinnest arm the study runs is a
# 50-row holdout block (n = 250 pooled, n = 500 unpooled), ~25 rows per
# treatment arm, so falling under 10 is a ~5-SD draw: if it happens, fail
# the run loudly rather than score with a near-empty fit.
MIN_ARM_ROWS <- 10L

#' The training rows of one treatment arm, with the MIN_ARM_ROWS guard
arm_rows <- function(W, filter_indices, w) {
  rows <- filter_indices & W == w
  if (sum(rows) < MIN_ARM_ROWS) {
    stop(sprintf(
      "only %d training rows with W = %d (MIN_ARM_ROWS = %d) - too few for a per-arm outcome model",
      sum(rows), w, MIN_ARM_ROWS
    ), call. = FALSE)
  }
  rows
}

prepare_xgb_matrices <- function(X, Y, W, filter_indices, training = TRUE) {
  if (training) {
    rows0 <- arm_rows(W, filter_indices, 0)
    rows1 <- arm_rows(W, filter_indices, 1)
    list(
      # per-arm outcome models - see the header
      X_Y0 = xgb.DMatrix(data = as.matrix(X[rows0, ]), label = Y[rows0]),
      X_Y1 = xgb.DMatrix(data = as.matrix(X[rows1, ]), label = Y[rows1]),
      X_W = xgb.DMatrix(
        data = as.matrix(X[filter_indices, ]),
        label = W[filter_indices]
      )
    )
  } else {
    # For prediction matrices
    list(
      X = xgb.DMatrix(data = as.matrix(X[filter_indices, ]))
    )
  }
}

fit_xgb <- function(grid, data_matrix) {
  cv_results <- run_xgb_cv(grid, data_matrix)
  trained_result <- train_best_xgb(grid, cv_results, data_matrix)
  return(trained_result$model)
}

run_xgb_cv <- function(grid, train_data) {
  cv_results <- apply(grid, 1, function(params) {
    xgb.cv(
      params = as.list(params),
      data = train_data,
      nrounds = 100,
      nfold = 5,
      early_stopping_rounds = 10,
      verbose = 0
    )
  })

  cv_results <- lapply(cv_results, function(params) {
    error_col <- if (grid$objective[1] == "reg:squarederror") {
      "test_rmse_mean"
    } else {
      "test_logloss_mean"
    }
    log_col <- params$evaluation_log[[error_col]]

    # best_iteration is top-level in xgboost <= 2 (the cluster) and under
    # $early_stop in xgboost >= 3 (local), so read both. which.min() of the log
    # is an exact equivalent: early stopping picks the minimum of this column.
    best_it <- params$early_stop$best_iteration   # xgboost >= 3
    if (length(best_it) != 1L || is.na(best_it)) {
      best_it <- params$best_iteration            # xgboost <= 2
    }
    if (length(best_it) != 1L || is.na(best_it)) {
      # neither location - fall back to the log itself. Only genuinely
      # degenerate if the log is empty too, which is what a training block
      # too small for a 5-fold inner CV would produce.
      if (length(log_col)) {
        best_it <- which.min(log_col)
      } else {
        best_it <- 1L
        warning(
          "xgb.cv() produced an empty evaluation log - training block too ",
          "small for a 5-fold inner CV. Falling back to nrounds = 1.",
          call. = FALSE
        )
      }
    }

    data.frame(
      nrounds = best_it,
      error = if (length(log_col)) min(log_col) else Inf
    )
  })

  return(do.call(rbind, cv_results))
}

train_best_xgb <- function(grid, cv_results, train_data) {
  best_idx <- which.min(cv_results$error)
  best_params <- as.list(grid[best_idx, ])
  best_nrounds <- cv_results$nrounds[best_idx]

  best_model <- xgb.train(
    params = best_params,
    data = train_data,
    nrounds = best_nrounds,
    verbose = 0
  )

  return(list(model = best_model, params = best_params, nrounds = best_nrounds))
}

fit_nuisance_xgb <- function(train_matrices, param_grids) {
  list(
    mu0 = fit_xgb(param_grids$continuous, train_matrices$X_Y0),
    mu1 = fit_xgb(param_grids$continuous, train_matrices$X_Y1),
    pi = fit_xgb(param_grids$binary, train_matrices$X_W)
  )
}

predict_nuisance_xgb <- function(models, pred_matrices) {
  data.frame(
    mu0_DR = predict(models$mu0, newdata = pred_matrices$X),
    mu1_DR = predict(models$mu1, newdata = pred_matrices$X),
    pi = predict(models$pi, newdata = pred_matrices$X)
  )
}

#' Fit-and-predict-on-the-same-block XGBoost nuisance estimation
#'
#' The `holdout` arm. The blocks are the candidate's own folds, from holdout_blocks(), not a
#' fresh draw: the point of the arm is that the nuisance sees only data the
#' candidate never saw. What it gives up in exchange is row-level honesty -
#' every row's own Y is in the model that predicts it, so mu_DR interpolates
#' it and phi's AIPW correction (Y - mu_DR) collapses toward zero. That is the
#' arm's defining property, not a defect to patch: cv_shared is the row-honest
#' arm, and the pair is what makes the cost of each visible.
run_xgb_holdout <- function(
  X,
  Y,
  W,
  fold_indices,
  fold_list,
  param_grids
) {
  nuisances <- lapply(fold_list, function(fold) {
    # the whole arm in one line: train and test are the same rows
    train_filter <- fold_indices %in% fold
    test_filter <- train_filter

    train_matrices <- prepare_xgb_matrices(
      X,
      Y,
      W,
      train_filter,
      training = TRUE
    )
    pred_matrices <- prepare_xgb_matrices(
      X,
      Y,
      W,
      test_filter,
      training = FALSE
    )

    models <- fit_nuisance_xgb(train_matrices, param_grids)

    predictions <- predict_nuisance_xgb(models, pred_matrices)
    predictions$row_index <- which(test_filter)

    return(predictions)
  })

  result <- unlist_order(nuisances)
  result <- calculate_pseudos(result, Y, W)

  return(result)
}

#' Leave-one-fold-out XGBoost nuisance estimation
#'
#' Serves the `cv_shared` arm, with the candidates' own fold_info handed over
#' by the caller.
run_xgb_cross_validation <- function(
  X,
  Y,
  W,
  fold_indices,
  fold_list,
  param_grids
) {
  nuisances <- lapply(fold_list, function(fold) {
    train_filter <- !(fold_indices %in% fold)
    test_filter <- !train_filter

    # Prepare data matrices
    train_matrices <- prepare_xgb_matrices(
      X,
      Y,
      W,
      train_filter,
      training = TRUE
    )
    pred_matrices <- prepare_xgb_matrices(
      X,
      Y,
      W,
      test_filter,
      training = FALSE
    )

    # Fit models
    models <- fit_nuisance_xgb(train_matrices, param_grids)

    # Make predictions
    predictions <- predict_nuisance_xgb(models, pred_matrices)
    predictions$row_index <- which(test_filter)

    return(predictions)
  })

  # Combine results and calculate pseudo-outcomes
  result <- unlist_order(nuisances)
  result <- calculate_pseudos(result, Y, W)

  return(result)
}

run_xgb_whole_dataset <- function(X, Y, W, param_grids) {
  n_obs <- nrow(X)
  all_indices <- rep(TRUE, n_obs)

  # Prepare data matrices
  train_matrices <- prepare_xgb_matrices(X, Y, W, all_indices, training = TRUE)
  pred_matrices <- prepare_xgb_matrices(X, Y, W, all_indices, training = FALSE)

  # Fit models
  models <- fit_nuisance_xgb(train_matrices, param_grids)

  # Make predictions
  result <- predict_nuisance_xgb(models, pred_matrices)
  result <- calculate_pseudos(result, Y, W)

  return(result)
}

#' Every XGBoost nuisance arm the caller asks for
#'
#' @param arms named list of arm specs - see nuisance_arm_spec(). Names become
#'   the names of the returned list, and therefore the arm half of
#'   me_metrics.R's score column names.
run_all_xgb_nuisance <- function(X, Y, W, arms, n_cores = 1) {
  param_grids <- create_param_grid(n_cores)

  lapply(arms, function(spec) {
    switch(
      spec$type,
      whole = run_xgb_whole_dataset(X, Y, W, param_grids),
      cv = run_xgb_cross_validation(
        X, Y, W, spec$fold_indices, spec$fold_list, param_grids
      ),
      holdout = run_xgb_holdout(
        X, Y, W, spec$fold_indices, spec$fold_list, param_grids
      ),
      stop("unknown nuisance arm type: ", spec$type)
    )
  })
}

# AutoML
# PBS Pro array jobs give $TMPDIR a path containing the literal
# "[<array index>]" (e.g. /var/tmp/pbs.<jobid>[<idx>].pbs-7/...) - a plain
# interactive session's $TMPDIR never has this. ice_root defaults to
# tempdir() (built from $TMPDIR), and H2O's Java-side path validator
# rejects "[" / "]" as illegal characters, so h2o.init() fails array-job-only
# with "Invalid ice_root ... Illegal character in path". Both setup and
# shutdown derive the same PID-scoped path, kept out of tempdir()/$TMPDIR
# entirely so it's immune to this (and to any other special characters PBS
# might put in a scratch path).
h2o_ice_root <- function() {
  file.path(
    "/var/tmp",
    sprintf("h2o_ice_%s_%d", Sys.getenv("USER", "u"), Sys.getpid())
  )
}

setup_h2o_cluster <- function(n_cores, mem, model_seed = NULL) {
  # Check if cluster is already running
  tryCatch(
    {
      h2o.clusterIsUp()
      cat("H2O cluster is already running\n")
    },
    error = function(e) {
      # Initialize H2O
      ice_root <- h2o_ice_root()
      dir.create(ice_root, recursive = TRUE, showWarnings = FALSE)
      h2o.init(
        port = parallelly::freePort(),
        nthreads = n_cores,
        max_mem_size = mem,
        min_mem_size = "1G",
        ice_root = ice_root
      )
    }
  )
}

# Set up shutdown handler
h2o_shutdown_check <- function() {
  cat("Ensuring H2O cluster shutdown...\n")
  tryCatch(
    {
      h2o.shutdown(prompt = FALSE)
      cat("H2O cluster shutdown successfully\n")
    },
    error = function(e) {
      cat("H2O cluster was already shutdown\n")
    }
  )
  unlink(h2o_ice_root(), recursive = TRUE, force = TRUE)
}

prepare_h2o_data <- function(X, Y, W, filter_indices, training = TRUE) {
  # Prepare different data combinations
  if (training) {
    YX <- cbind(Y, X)
    WX <- cbind(W = as.factor(W), X)

    list(
      # per-arm outcome models - see the header
      YX0 = as.h2o(YX[arm_rows(W, filter_indices, 0), ]),
      YX1 = as.h2o(YX[arm_rows(W, filter_indices, 1), ]),
      WX = as.h2o(WX[filter_indices, ])
    )
  } else {
    list(
      X = as.h2o(X[filter_indices, ])
    )
  }
}

fit_automl_models <- function(
  train_data,
  model_seed,
  max_models = 20,
  exclude_algos = c("DeepLearning", "XGBoost")
) {
  list(
    mu0 = h2o.automl(
      y = "Y",
      training_frame = train_data$YX0,
      nfolds = 0,
      max_models = max_models,
      exclude_algos = exclude_algos,
      seed = model_seed
    ),
    mu1 = h2o.automl(
      y = "Y",
      training_frame = train_data$YX1,
      nfolds = 0,
      max_models = max_models,
      exclude_algos = exclude_algos,
      seed = model_seed
    ),
    pi = h2o.automl(
      y = "W",
      training_frame = train_data$WX,
      nfolds = 0,
      max_models = max_models,
      exclude_algos = exclude_algos,
      seed = model_seed
    )
  )
}

predict_automl_models <- function(models, pred_data) {
  pi_pred <- h2o.predict(models$pi, newdata = pred_data$X) %>% as.data.frame()

  data.frame(
    mu0_DR = h2o.predict(models$mu0, newdata = pred_data$X) %>% as.vector(),
    mu1_DR = h2o.predict(models$mu1, newdata = pred_data$X) %>% as.vector(),
    pi = pi_pred$p1
  )
}

#' Fit-and-predict-on-the-same-block H2O AutoML nuisance estimation
#'
#' The AutoML half of the `holdout` arm - see run_xgb_holdout() for what the
#' arm is for and what it trades away.
#'
#' W goes to calculate_pseudos() unchanged, as in the other arms. Don't "undo"
#' H2O's factor with as.numeric(W) - 1: prepare_h2o_data() factors only its own
#' copy, so W here is still 0/1 and that would make it -1/0.
run_automl_holdout <- function(
  X,
  Y,
  W,
  fold_indices,
  fold_list,
  model_seed,
  max_models = 20,
  exclude_algos = c("DeepLearning", "XGBoost")
) {
  nuisances <- lapply(fold_list, function(fold) {
    train_filter <- fold_indices %in% fold
    test_filter <- train_filter

    train_data <- prepare_h2o_data(X, Y, W, train_filter, training = TRUE)
    pred_data <- prepare_h2o_data(X, Y, W, test_filter, training = FALSE)

    models <- fit_automl_models(
      train_data,
      model_seed,
      max_models,
      exclude_algos
    )

    predictions <- predict_automl_models(models, pred_data)
    predictions$row_index <- which(test_filter)

    h2o.removeAll()

    return(predictions)
  })

  result <- unlist_order(nuisances)
  result <- calculate_pseudos(result, Y, W)

  return(result)
}

#' Leave-one-fold-out H2O AutoML nuisance estimation
#'
#' Serves the `cv_shared` arm, as the XGBoost version does.
run_automl_cross_validation <- function(
  X,
  Y,
  W,
  fold_indices,
  fold_list,
  model_seed,
  max_models = 20,
  exclude_algos = c("DeepLearning", "XGBoost")
) {
  nuisances <- lapply(fold_list, function(fold) {
    train_filter <- !(fold_indices %in% fold)
    test_filter <- !train_filter

    # Prepare data
    train_data <- prepare_h2o_data(X, Y, W, train_filter, training = TRUE)
    pred_data <- prepare_h2o_data(X, Y, W, test_filter, training = FALSE)

    # Fit models
    models <- fit_automl_models(
      train_data,
      model_seed,
      max_models,
      exclude_algos
    )

    # Make predictions
    predictions <- predict_automl_models(models, pred_data)
    predictions$row_index <- which(test_filter)

    # Clean up H2O objects for this fold
    h2o.removeAll()

    return(predictions)
  })

  # Combine results and calculate pseudo-outcomes
  result <- unlist_order(nuisances)
  result <- calculate_pseudos(result, Y, W)

  return(result)
}

run_automl_whole_dataset <- function(
  X,
  Y,
  W,
  model_seed,
  max_models = 20,
  exclude_algos = c("DeepLearning", "XGBoost")
) {
  n_obs <- nrow(X)
  all_indices <- rep(TRUE, n_obs)

  train_data <- prepare_h2o_data(X, Y, W, all_indices, training = TRUE)
  pred_data <- prepare_h2o_data(X, Y, W, all_indices, training = FALSE)

  # Fit models
  models <- fit_automl_models(train_data, model_seed, max_models, exclude_algos)

  # Make predictions
  result <- predict_automl_models(models, pred_data)
  result <- calculate_pseudos(result, Y, W)

  return(result)
}

#' Every AutoML nuisance arm the caller asks for
#'
#' The cluster is started once for all arms and shut down once at the end -
#' not per arm. Each h2o.init() is a JVM launch, and the arm list is now three
#' or four long rather than two.
#'
#' @param arms named list of arm specs - see nuisance_arm_spec()
run_all_automl_nuisance <- function(
  X,
  Y,
  W,
  arms,
  n_cores = 1,
  mem = "5G",
  model_seed,
  max_models = 20,
  exclude_algos = c("DeepLearning", "XGBoost")
) {
  # Setup H2O cluster
  setup_h2o_cluster(n_cores = n_cores, mem = mem, model_seed = model_seed)
  on.exit(h2o_shutdown_check(), add = TRUE)

  lapply(arms, function(spec) {
    switch(
      spec$type,
      whole = run_automl_whole_dataset(
        X, Y, W, model_seed, max_models, exclude_algos
      ),
      cv = run_automl_cross_validation(
        X, Y, W, spec$fold_indices, spec$fold_list,
        model_seed, max_models, exclude_algos
      ),
      holdout = run_automl_holdout(
        X, Y, W, spec$fold_indices, spec$fold_list,
        model_seed, max_models, exclude_algos
      ),
      stop("unknown nuisance arm type: ", spec$type)
    )
  })
}

#' One arm's specification
#'
#' @param type "whole", "cv" (leave-one-fold-out) or "holdout" (fit and
#'   predict on the same block)
#' @param folds for cv/holdout, a list of fold_indices and fold_list as
#'   returned by split_folds() or holdout_blocks()
nuisance_arm_spec <- function(type, folds = NULL) {
  type <- match.arg(type, c("whole", "cv", "holdout"))
  if (type == "whole") return(list(type = type))
  stopifnot(!is.null(folds), !is.null(folds$fold_indices), !is.null(folds$fold_list))
  list(type = type, fold_indices = folds$fold_indices, fold_list = folds$fold_list)
}

#' Both nuisance-evaluation pipelines, over an arbitrary set of arms
#'
#' The one entry point, used by me_analysis.R, me_strategies.R and me_split.R.
#'
#' @param arms named list of nuisance_arm_spec()s. Names become the arm half
#'   of me_metrics.R's score column names, so they should come from
#'   me_config.R's NUISANCE_ARMS rather than being typed at the call site.
#' @param pipelines which estimators to run. Restricting to "xgb" is what lets
#'   me_testing.R exercise the arm plumbing locally without a Java runtime.
run_nuisance_arms <- function(
  X,
  Y,
  W,
  arms,
  n_cores,
  mem,
  model_seed,
  pipelines = c("xgb", "automl")
) {
  out <- list()

  if ("xgb" %in% pipelines) {
    out$xgb <- run_all_xgb_nuisance(X, Y, W, arms, n_cores = n_cores)
  }
  if ("automl" %in% pipelines) {
    out$automl <- run_all_automl_nuisance(
      X, Y, W, arms,
      n_cores = n_cores, mem = mem, model_seed = model_seed
    )
  }

  out
}
