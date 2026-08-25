##########
# title: nuisance-evaluation pipelines - model evaluation study
##########
# A second, INDEPENDENT nuisance-estimation pipeline used only to score the 9
# candidate models from me_models.R against each other and against true PEHE
# - not to fit them. No other study in this repo needs this concept, so it
# doesn't map onto config/dgms/models/analysis/check/collect/metrics and gets
# its own file, the same way crossfitting/ has files beyond that 7-role
# floor (cf_testing.R, cf_profile.R, ...)
#
# Two independent estimators of the same 4 nuisance targets (mu0_T, mu1_T,
# mu_DR, pi): XGBoost with a hand-tuned CV grid, and H2O AutoML. Each runs
# under a set of ARMS differing in one thing only: what data the nuisance sees
# relative to the data the candidate it is scoring was trained on. The arm
# names and what they mean live in me_config.R's NUISANCE_ARMS.
#
# WHAT CHANGED, AND WHY THE OLD RATIONALE IS GONE. This file used to document
# the fold_indices passed in as "an INDEPENDENT draw from the one used to fit
# the 9 candidates", on the grounds that sharing one split would fit a
# candidate's tau_hat and this pipeline's nuisance at the same row on the
# identical training set, correlating their errors. That claim does not
# survive the arithmetic:
#
#   - Row-level honesty - whether row i's own (Y_i, W_i) entered the model
#     predicting at row i - holds under BOTH a shared and an independent draw.
#     Neither lets a nuisance see the row it predicts. The second draw buys
#     nothing on the axis it was justified by.
#   - What it does change is the overlap of the two training sets for row i,
#     from 100% down to (V-1)/V = 90%. It removes a tenth of the dependence.
#
# crossfitting/README.md makes the general version of this point ("Re-
# randomising the stage-2 split cannot remove that dependence"), which is why
# scf_scf_new is an arm to be TESTED there rather than the default. So the
# independent draw is retired as a default here, kept as the `cv_indep` arm
# because it is already computed, and set beside `cv_shared` - which uses the
# candidates' own folds - so the comparison is empirical rather than asserted.
#
# The `holdout` arm restores the `infold` (resubstitution) branch that was
# removed from this file: fit AND predict on the fold the candidate held out.
# It is deliberately not row-level honest - it is the arm that is fully
# decoupled from the candidate instead, and the two properties cannot both be
# had from a single small block. See the comments above the two holdout
# functions for what the removal took out and which part of it must NOT come
# back.
require(h2o)
require(xgboost)

unlist_order <- function(reslist) {
  df <- do.call(rbind, reslist)
  df <- df %>%
    arrange(row_index) %>%
    select(-row_index)
  return(df)
}

calculate_pseudos <- function(df, Y, W) {
  df <- df %>%
    mutate(
      mu_DR = mu1_DR * W + (1 - W) * mu0_DR
    )
  df <- df %>%
    mutate(
      tau_T = mu1_T - mu0_T,
      phi = mu1_DR - mu0_DR + ((Y - mu_DR) * (W - pi)) / (pi * (1 - pi)),
      # phi05 (propensity fixed at 0.5) is currently unused by any score in
      # me_metrics.R - kept since it's one extra column from data already in
      # scope, not worth deleting the way a whole unused function is.
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

prepare_xgb_matrices <- function(X, Y, W, filter_indices, training = TRUE) {
  if (training) {
    list(
      X_Y = xgb.DMatrix(
        data = as.matrix(X[filter_indices, ]),
        label = Y[filter_indices]
      ),
      X_0_Y = xgb.DMatrix(
        data = as.matrix(X[filter_indices & W == 0, ]),
        label = Y[filter_indices & W == 0]
      ),
      X_1_Y = xgb.DMatrix(
        data = as.matrix(X[filter_indices & W == 1, ]),
        label = Y[filter_indices & W == 1]
      ),
      XW_Y = xgb.DMatrix(
        data = as.matrix(cbind(W = W[filter_indices], X[filter_indices, ])),
        label = Y[filter_indices]
      ),
      X_W = xgb.DMatrix(
        data = as.matrix(X[filter_indices, ]),
        label = W[filter_indices]
      )
    )
  } else {
    # For prediction matrices
    list(
      X = xgb.DMatrix(data = as.matrix(X[filter_indices, ])),
      X0 = xgb.DMatrix(data = as.matrix(cbind(W = 0, X[filter_indices, ]))),
      X1 = xgb.DMatrix(data = as.matrix(cbind(W = 1, X[filter_indices, ])))
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

    # WHERE best_iteration LIVES DEPENDS ON THE xgboost VERSION. Up to
    # xgboost 2.x it is a top-level element of the xgb.cv() result, which is
    # what this code originally read. xgboost 3.x moved it into an
    # `early_stop` sub-list, so `params$best_iteration` is integer(0) there -
    # for EVERY fit, at every data size, not just small ones. That made
    # data.frame() error with "arguments imply differing number of rows: 0, 1"
    # and would take a whole array task down.
    #
    # Local R has xgboost 3.2.1.1; the cluster's R/4.3.2 module has an older
    # one (the 358 completed runs prove it, since they went through this exact
    # line). So this is the same class of local/cluster version mismatch as
    # SL2's - see README.md - and the fix has to work under both.
    #
    # The which.min() fallback is a genuine equivalent, not a degradation:
    # early stopping selects the minimum of the very column being scanned
    # here, and the two agree exactly where both are available. So this
    # changes no number the cluster already produced.
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
    mu0_T = fit_xgb(param_grids$continuous, train_matrices$X_0_Y),
    mu1_T = fit_xgb(param_grids$continuous, train_matrices$X_1_Y),
    mu_DR = fit_xgb(param_grids$continuous, train_matrices$XW_Y),
    pi = fit_xgb(param_grids$binary, train_matrices$X_W)
  )
}

predict_nuisance_xgb <- function(models, pred_matrices) {
  data.frame(
    mu0_T = predict(models$mu0_T, newdata = pred_matrices$X),
    mu1_T = predict(models$mu1_T, newdata = pred_matrices$X),
    mu0_DR = predict(models$mu_DR, newdata = pred_matrices$X0),
    mu1_DR = predict(models$mu_DR, newdata = pred_matrices$X1),
    pi = predict(models$pi, newdata = pred_matrices$X)
  )
}

#' Fit-and-predict-on-the-same-block XGBoost nuisance estimation
#'
#' The `holdout` arm. This is the `infold` branch that used to live inside
#' run_xgb_cross_validation() behind a `method` argument and was removed; the
#' note left behind said a future re-add "can restore the branch above
#' unchanged", and that is what this is. XGBoost's version had no companion
#' bug (AutoML's did - see run_automl_holdout()).
#'
#' The blocks are the candidate's OWN folds, from holdout_blocks(), not a
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
#' Serves both `cv_indep` and `cv_shared` - the two differ only in which fold
#' vector the caller hands over (an independent draw vs. the candidate's own
#' fold_info), not in anything this function does.
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
#' Fixes a bug found while porting: the old sim_eval.R called this without
#' n_cores, so XGBoost's own nthread grid parameter silently defaulted to 1
#' regardless of the n_cores the script set elsewhere. run_all_nuisance_pipelines()
#' below now always passes it through explicitly.
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
    YWX <- cbind(Y, W, X)
    WX <- cbind(W = as.factor(W), X)

    list(
      YX_0 = as.h2o(YX[filter_indices & W == 0, ]),
      YX_1 = as.h2o(YX[filter_indices & W == 1, ]),
      YWX = as.h2o(YWX[filter_indices, ]),
      WX = as.h2o(WX[filter_indices, ])
    )
  } else {
    list(
      X = as.h2o(X[filter_indices, ]),
      X0 = as.h2o(cbind(W = 0, X)[filter_indices, ]),
      X1 = as.h2o(cbind(W = 1, X)[filter_indices, ])
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
    mu0_T = h2o.automl(
      y = "Y",
      training_frame = train_data$YX_0,
      nfolds = 0,
      max_models = max_models,
      exclude_algos = exclude_algos,
      seed = model_seed
    ),
    mu1_T = h2o.automl(
      y = "Y",
      training_frame = train_data$YX_1,
      nfolds = 0,
      max_models = max_models,
      exclude_algos = exclude_algos,
      seed = model_seed
    ),
    mu_DR = h2o.automl(
      y = "Y",
      training_frame = train_data$YWX,
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
    mu0_T = h2o.predict(models$mu0_T, newdata = pred_data$X) %>% as.vector(),
    mu1_T = h2o.predict(models$mu1_T, newdata = pred_data$X) %>% as.vector(),
    mu0_DR = h2o.predict(models$mu_DR, newdata = pred_data$X0) %>% as.vector(),
    mu1_DR = h2o.predict(models$mu_DR, newdata = pred_data$X1) %>% as.vector(),
    pi = pi_pred$p1
  )
}

#' Fit-and-predict-on-the-same-block H2O AutoML nuisance estimation
#'
#' The AutoML half of the `holdout` arm - see run_xgb_holdout() for what the
#' arm is for and what it trades away.
#'
#' THE PART OF THE OLD `infold` BRANCH THAT IS DELIBERATELY NOT RESTORED. The
#' removed code was not just the train/test filter below; it shipped with a
#' W-recoding line:
#'
#'   W_final <- if (method == "infold") as.numeric(W) - 1 else W
#'   result <- calculate_pseudos(result, Y, W_final)
#'
#' prepare_h2o_data()'s `WX <- cbind(W = as.factor(W), X)` turns W into a
#' factor for H2O's classification target, and as.numeric() on a 2-level
#' factor gives 1/2 rather than 0/1 - so the line reads as an attempt to undo
#' that. But `W` in THIS function's scope is the caller's original numeric
#' vector; it is never itself converted to a factor. as.numeric(W) - 1 would
#' therefore take a plain 0/1 vector to -1/0, and calculate_pseudos()'s mu_DR
#' and phi both assume W in {0,1}. The old note flagged this as unverified
#' rather than fixed, because the pipeline had never completed a run; it has
#' now, and the reading holds. The filter comes back, the recoding does not -
#' W is passed through to calculate_pseudos() unchanged, exactly as the cv and
#' whole arms already do.
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
#' Serves both `cv_indep` and `cv_shared`, as the XGBoost version does - only
#' the fold vector handed in differs.
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
#' The cluster is started once for ALL arms and shut down once at the end -
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
  # was on.exit(h2o_shutdown_check, add = TRUE) - missing parens meant the
  # cleanup handler was registered but never actually called on exit, so the
  # "always shut the cluster down, even on error" safety net never fired.
  # The explicit h2o.shutdown() this function used to end with on its happy
  # path is now redundant (h2o_shutdown_check()'s own tryCatch makes it safe
  # to call twice, but there's no need to) and has been removed.
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
#' The general entry point, used by me_strategies.R and me_split.R.
#' run_all_nuisance_pipelines() below is the legacy two-arm wrapper.
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

#' The original two-arm pipeline, unchanged in behaviour
#'
#' Kept at its old signature ON PURPOSE. me_analysis.R has already produced
#' 358 of 360 runs against it, and rerunning one of the two failures must
#' still write the same `cv`/`whole` shape those files have - me_strategies.R
#' is what renames `cv` to `cv_indep` and adds the rest, and it has to find
#' the same thing in every input file regardless of when that file was
#' written.
run_all_nuisance_pipelines <- function(
  X,
  Y,
  W,
  fold_indices,
  fold_list,
  n_cores,
  mem,
  model_seed
) {
  folds <- list(fold_indices = fold_indices, fold_list = fold_list)

  run_nuisance_arms(
    X, Y, W,
    arms = list(
      cv = nuisance_arm_spec("cv", folds),
      whole = nuisance_arm_spec("whole")
    ),
    n_cores = n_cores,
    mem = mem,
    model_seed = model_seed
  )
}
