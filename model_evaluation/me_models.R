##########
# title: candidate CATE-learner fitting - model evaluation study
##########
# The 9 candidate configurations this study compares: 3 random-forest
# hyperparameter sets, 3 elastic-net (glmnet) sets, 3 SuperLearner library
# sets. Each is fit as its own double-crossfit DR-learner - out of scope to
# change, this is what the study evaluates.
#
# This is a *different* estimator family from R/cate_models.R's DR-learner
# (e.g. fit_glmnet wraps a single learner via create.Learner("SL.glmnet", ...)
# inside SuperLearner(), a code path R/cate_models.R doesn't have) - not
# unified into R/ here, since that would be a redesign, not a plumbing port.
#
# collate_predictions() is sourced from R/utils.R rather than duplicated
# locally (the old model_utils.R's copy was logically identical, just
# %>%-piped instead of nested).

require(future.apply)
require(ranger)
require(SuperLearner)

# hyper parameter specifications for CATE learners
create_rf_hyperparams <- function(
  mtry = NULL,
  max.depth = NULL,
  sample.fraction = 0.5,
  replace = FALSE
) {
  list(
    mtry = mtry,
    max.depth = max.depth,
    sample.fraction = sample.fraction,
    replace = replace
  )
}

create_net_hyperparams <- function(alpha) {
  list(alpha = alpha)
}

create_SL_hyperparams <- function(method, SL.library) {
  list(method = method, SL.library = SL.library)
}

# general fitting functions
collect_tau <- function(tau_list, fold_list, fold_indices) {
  tau <- rep(NA, length(fold_indices))
  for (fold in fold_list) {
    tau[fold_indices == fold] <- tau_list[[fold]]
  }
  return(tau)
}

# fitting functions for the random forest
fit_rf <- function(
  Y,
  W,
  X,
  hyper_list = list(),
  fold_indices,
  fold_list,
  fold_pairs
) {
  hyper_list$write.forest <- TRUE

  nuisances <- future_lapply(
    fold_pairs,
    crossfit_double_RF,
    Y,
    W,
    X,
    hyper_list,
    fold_indices,
    future.seed = TRUE
  )

  predictions <- c("po", "Y0.hat", "Y1.hat", "W.hat")
  matrices <- setNames(
    lapply(predictions, function(pred) {
      collate_predictions(fold_list, fold_pairs, fold_indices, nuisances, pred)
    }),
    predictions
  )

  tau <- fit_tau_rf(X, fold_list, fold_indices, matrices$po, hyper_list)

  c(list(tau = tau), matrices)
}

crossfit_double_RF <- function(fold_pair, Y, W, X, hyper_list, fold_indices) {
  train_filter <- !(fold_indices %in% fold_pair)
  test_filter <- !train_filter

  fit_model <- function(formula, data) {
    hyper <- hyper_list
    hyper$data <- data
    hyper$formula <- formula
    do.call(ranger, hyper)
  }

  Y_model <- fit_model(
    Y ~ .,
    data.frame(Y = Y[train_filter], W = W[train_filter], X[train_filter, ])
  )
  W_model <- fit_model(
    W ~ .,
    data.frame(W = W[train_filter], X[train_filter, ])
  )

  test_data <- data.frame(W, X)[test_filter, ]
  Y0.hat <- predict(Y_model, mutate(test_data, W = 0L))$predictions
  Y1.hat <- predict(Y_model, mutate(test_data, W = 1L))$predictions
  W.hat <- predict(W_model, X[test_filter, ])$predictions

  W_test <- W[test_filter]
  Y_test <- Y[test_filter]
  Y.hat <- W_test * Y1.hat + (1 - W_test) * Y0.hat

  po <- Y1.hat -
    Y0.hat +
    ((Y_test - Y.hat) * (W_test - W.hat)) / (W.hat * (1 - W.hat))

  list(po = po, Y0.hat = Y0.hat, Y1.hat = Y1.hat, W.hat = W.hat)
}

fit_tau_rf <- function(X, fold_list, fold_indices, po_matrix, hyper_list) {
  tau <- rep(NA, length(fold_indices))

  tau_list <- future_lapply(
    fold_list,
    function(fold) {
      train_filter <- !(fold_indices %in% fold)
      test_filter <- !train_filter

      hyper_po <- hyper_list
      hyper_po$data <- data.frame(
        po = po_matrix[train_filter, fold],
        X[train_filter, ]
      )
      hyper_po$formula <- po ~ .

      po_model <- do.call(ranger, hyper_po)
      predict(po_model, X[test_filter, ])$predictions
    },
    future.seed = TRUE
  )

  collect_tau(tau_list, fold_list, fold_indices)
}

# fitting functions for GLMnet
fit_glmnet <- function(
  Y,
  W,
  X,
  hyper_list = list(),
  fold_indices,
  fold_list,
  fold_pairs
) {
  nuisances <- future_lapply(
    fold_pairs,
    crossfit_double_glmnet,
    Y,
    W,
    X,
    hyper_list,
    fold_indices,
    future.seed = TRUE
  )

  predictions <- c("po", "Y0.hat", "Y1.hat", "W.hat")
  matrices <- setNames(
    lapply(predictions, function(pred) {
      collate_predictions(fold_list, fold_pairs, fold_indices, nuisances, pred)
    }),
    predictions
  )

  tau <- fit_tau_glmnet(X, fold_list, fold_indices, matrices$po, hyper_list)

  c(list(tau = tau), matrices)
}

crossfit_double_glmnet <- function(
  fold_pair,
  Y,
  W,
  X,
  hyper_list,
  fold_indices
) {
  custom_net <- create.Learner("SL.glmnet", params = hyper_list)

  train_filter <- !(fold_indices %in% fold_pair)
  test_filter <- !train_filter

  fit_model <- function(y, x, family) {
    SuperLearner(
      y,
      x,
      family = family,
      SL.library = custom_net$names,
      method = "method.CC_LS"
    )
  }

  Y_model <- fit_model(Y[train_filter], cbind(W, X)[train_filter, ], gaussian())
  W_model <- fit_model(W[train_filter], X[train_filter, ], binomial())

  test_data <- data.frame(W, X)[test_filter, ]
  Y0.hat <- predict(Y_model, mutate(test_data, W = 0L))$pred
  Y1.hat <- predict(Y_model, mutate(test_data, W = 1L))$pred
  W.hat <- predict(W_model, X[test_filter, ])$pred

  W_test <- W[test_filter]
  Y_test <- Y[test_filter]
  Y.hat <- W_test * Y1.hat + (1 - W_test) * Y0.hat

  po <- Y1.hat -
    Y0.hat +
    ((Y_test - Y.hat) * (W_test - W.hat)) / (W.hat * (1 - W.hat))

  list(po = po, Y0.hat = Y0.hat, Y1.hat = Y1.hat, W.hat = W.hat)
}

fit_tau_glmnet <- function(X, fold_list, fold_indices, po_matrix, hyper_list) {
  custom_net <- create.Learner("SL.glmnet", params = hyper_list)

  tau_list <- future_lapply(
    fold_list,
    function(fold) {
      train_filter <- !(fold_indices %in% fold)
      test_filter <- !train_filter

      po_model <- SuperLearner(
        po_matrix[train_filter, fold],
        X[train_filter, ],
        family = gaussian(),
        SL.library = custom_net$names,
        method = "method.CC_LS"
      )

      predict(po_model, X[test_filter, ])$pred
    },
    future.seed = TRUE
  )

  collect_tau(tau_list, fold_list, fold_indices)
}

# fitting functions for SuperLearner
fit_SL <- function(
  Y,
  W,
  X,
  hyper_list = list(),
  fold_indices,
  fold_list,
  fold_pairs
) {
  nuisances <- future_lapply(
    fold_pairs,
    crossfit_double_SL,
    Y,
    W,
    X,
    hyper_list,
    fold_indices,
    future.seed = TRUE
  )

  predictions <- c("po", "Y0.hat", "Y1.hat", "W.hat")
  matrices <- setNames(
    lapply(predictions, function(pred) {
      collate_predictions(fold_list, fold_pairs, fold_indices, nuisances, pred)
    }),
    predictions
  )

  tau <- fit_tau_SL(X, fold_list, fold_indices, matrices$po, hyper_list)

  c(list(tau = tau), matrices)
}

crossfit_double_SL <- function(fold_pair, Y, W, X, hyper_list, fold_indices) {
  train_filter <- !(fold_indices %in% fold_pair)
  test_filter <- !train_filter

  fit_model <- function(y, x, family) {
    hyper <- hyper_list
    hyper$family <- family
    hyper$Y <- y
    hyper$X <- x
    do.call(SuperLearner, hyper)
  }

  Y_model <- fit_model(Y[train_filter], cbind(W, X)[train_filter, ], gaussian())
  W_model <- fit_model(W[train_filter], X[train_filter, ], binomial())

  test_data <- data.frame(W, X)[test_filter, ]
  Y0.hat <- predict(Y_model, mutate(test_data, W = 0L))$pred
  Y1.hat <- predict(Y_model, mutate(test_data, W = 1L))$pred
  W.hat <- predict(W_model, X[test_filter, ])$pred

  W_test <- W[test_filter]
  Y_test <- Y[test_filter]
  Y.hat <- W_test * Y1.hat + (1 - W_test) * Y0.hat

  po <- Y1.hat -
    Y0.hat +
    ((Y_test - Y.hat) * (W_test - W.hat)) / (W.hat * (1 - W.hat))

  list(po = po, Y0.hat = Y0.hat, Y1.hat = Y1.hat, W.hat = W.hat)
}

fit_tau_SL <- function(X, fold_list, fold_indices, po_matrix, hyper_list) {
  tau_list <- future_lapply(
    fold_list,
    function(fold) {
      train_filter <- !(fold_indices %in% fold)
      test_filter <- !train_filter

      hyper_po <- hyper_list
      hyper_po$family <- gaussian()
      hyper_po$X <- X[train_filter, ]
      hyper_po$Y <- po_matrix[train_filter, fold]

      po_model <- do.call(SuperLearner, hyper_po)
      predict(po_model, X[test_filter, ])$pred
    },
    future.seed = TRUE
  )

  collect_tau(tau_list, fold_list, fold_indices)
}

#' Fit all 9 candidate CATE-learner configurations
#'
#' Replaces the old sim_eval.R's list2env()/.GlobalEnv/grep()/assign() dance
#' with an ordinary named list. Names here must match CANDIDATE_MODELS in
#' me_config.R exactly - me_metrics.R relies on it via
#' R/metrics.R::compute_metrics()'s intersect(names(sim_res), models).
#'
#' @return named list, one element per candidate, each
#'   list(tau=, po=, Y0.hat=, Y1.hat=, W.hat=)
run_all_candidate_models <- function(Y, W, X, fold_indices, fold_list, fold_pairs) {
  # rf2/rf3's mtry is scaled to ncol(X) rather than fixed at 30/10. Those
  # fixed values came from the old benchtm-based prototype, which had a much
  # wider covariate set than R/dgm_scenarios.R produces (7-9 columns across
  # this study's scenarios) - left as literals they make ranger error out
  # ("mtry can not be larger than number of variables in data"), which is
  # how this was actually found (me_testing.R full, check 4). rf2 uses every
  # covariate each split (the high-mtry extreme, well-defined at any p);
  # rf3 uses about half, paired with its original shallower max.depth = 3 -
  # keeping the three RF configs meaningfully separated (~sqrt(p) < p/2 < p)
  # regardless of how many covariates a given scenario produces.
  p <- ncol(X)

  hyperparams <- list(
    rf1 = list(), # ranger defaults (mtry = floor(sqrt(p)))
    rf2 = create_rf_hyperparams(mtry = p, max.depth = 5),
    rf3 = create_rf_hyperparams(mtry = max(2, ceiling(p / 2)), max.depth = 3),

    net1 = create_net_hyperparams(alpha = 1), # lasso
    net2 = create_net_hyperparams(alpha = 0), # ridge
    net3 = create_net_hyperparams(alpha = 0.5), # elastic net

    SL1 = create_SL_hyperparams(
      method = "method.CC_LS",
      SL.library = list("SL.glmnet", "SL.ranger", "SL.earth", "SL.gam", "SL.mean")
    ),
    SL2 = create_SL_hyperparams(
      method = "method.CC_LS",
      SL.library = list(
        "SL.glmnet", "SL.xgboost", "SL.cforest", "SL.earth", "SL.gam", "SL.mean"
      )
    ),
    SL3 = create_SL_hyperparams(
      method = "method.CC_LS",
      SL.library = list("SL.svm", "SL.nnet", "SL.mean")
    )
  )

  fitter_for <- function(nm) {
    if (grepl("^rf", nm)) fit_rf
    else if (grepl("^net", nm)) fit_glmnet
    else fit_SL
  }

  out <- lapply(names(hyperparams), function(nm) {
    fitter_for(nm)(Y, W, X, hyper_list = hyperparams[[nm]], fold_indices, fold_list, fold_pairs)
  })

  setNames(out, names(hyperparams))
}
