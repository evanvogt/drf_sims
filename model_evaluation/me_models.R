##########
# title: candidate CATE-learner fitting - model evaluation study
##########
# The 9 candidate configurations this study compares: 3 random-forest
# hyperparameter sets, 3 elastic-net (glmnet) sets, 3 SuperLearner library
# sets. Each is fit as its own single-crossfit ("scf_scf") DR-learner:
# leave-one-fold-out nuisances feeding a stage-2 regression crossfit over the
# SAME folds - not the double-crossfit-over-fold-pairs scheme this file used
# to implement.
#
# Brought in line with crossfitting/'s comparison of alternatives against
# double crossfitting (fitting nuisances over all C(V,2) fold pairs - 45 fits
# at V=10 rather than 10), which is why R/cate_models.R's production
# estimators moved off it: forests to whole-sample OOB, causal forest to
# grf's own internal cross-fitting, SuperLearner to scf_scf (nuisance_sl /
# stage_2_sl). See crossfitting/README.md and R/cate_models.R's header.
#
# Single crossfitting throughout here, not a per-learner OOB/scf split, even
# though rf1-3 (ranger) could take OOB predictions natively. net1-3 and SL1-3
# have no OOB analogue - glmnet has no bagging, and SuperLearner's internal CV
# selects ensemble weights rather than producing an honest prediction of a
# training row (crossfitting/README.md's "DR-learner, SuperLearner" section
# states this directly) - and this study's question ("do cheap proxy losses
# rank 9 candidate models the way true PEHE would?") needs all 9 candidates on
# one shared honesty regime, or the ranking confounds learner choice with
# crossfitting scheme.
#
# This is a *different* estimator family from R/cate_models.R's DR-learner
# (e.g. fit_glmnet wraps a single learner via create.Learner("SL.glmnet", ...)
# inside SuperLearner(), a code path R/cate_models.R doesn't have) - not
# unified into R/ here, since that would be a redesign, not a plumbing port.
#
# scatter_folds() is sourced from R/utils.R rather than duplicated locally -
# it is the single-crossfit reassembly helper R/cate_models.R::nuisance_sl
# uses, replacing this file's former use of collate_predictions() (the
# double-crossfit n x V matrix assembler, still used by crossfitting/).

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
  fold_list
) {
  hyper_list$write.forest <- TRUE

  nuisances <- future_lapply(
    fold_list,
    crossfit_single_RF,
    Y,
    W,
    X,
    hyper_list,
    fold_indices,
    future.seed = TRUE
  )

  matrices <- scatter_folds(
    nuisances, fold_indices, c("po", "Y0.hat", "Y1.hat", "W.hat")
  )

  tau <- fit_tau_rf(X, fold_list, fold_indices, matrices$po, hyper_list)

  c(list(tau = tau), matrices)
}

crossfit_single_RF <- function(fold, Y, W, X, hyper_list, fold_indices) {
  train_filter <- fold_indices != fold
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
  W.hat <- trim_ps(predict(W_model, X[test_filter, ])$predictions)

  W_test <- W[test_filter]
  Y_test <- Y[test_filter]
  Y.hat <- W_test * Y1.hat + (1 - W_test) * Y0.hat

  po <- Y1.hat -
    Y0.hat +
    ((Y_test - Y.hat) * (W_test - W.hat)) / (W.hat * (1 - W.hat))

  list(fold = fold, po = po, Y0.hat = Y0.hat, Y1.hat = Y1.hat, W.hat = W.hat)
}

fit_tau_rf <- function(X, fold_list, fold_indices, po, hyper_list) {
  tau_list <- future_lapply(
    fold_list,
    function(fold) {
      train_filter <- fold_indices != fold
      test_filter <- !train_filter

      hyper_po <- hyper_list
      hyper_po$data <- data.frame(
        po = po[train_filter],
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
  fold_list
) {
  nuisances <- future_lapply(
    fold_list,
    crossfit_single_glmnet,
    Y,
    W,
    X,
    hyper_list,
    fold_indices,
    future.seed = TRUE
  )

  matrices <- scatter_folds(
    nuisances, fold_indices, c("po", "Y0.hat", "Y1.hat", "W.hat")
  )

  tau <- fit_tau_glmnet(X, fold_list, fold_indices, matrices$po, hyper_list)

  c(list(tau = tau), matrices)
}

crossfit_single_glmnet <- function(
  fold,
  Y,
  W,
  X,
  hyper_list,
  fold_indices
) {
  custom_net <- create.Learner("SL.glmnet", params = hyper_list)

  train_filter <- fold_indices != fold
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
  Y0.hat <- as.numeric(predict(Y_model, mutate(test_data, W = 0L))$pred)
  Y1.hat <- as.numeric(predict(Y_model, mutate(test_data, W = 1L))$pred)
  W.hat <- trim_ps(as.numeric(predict(W_model, X[test_filter, ])$pred))

  W_test <- W[test_filter]
  Y_test <- Y[test_filter]
  Y.hat <- W_test * Y1.hat + (1 - W_test) * Y0.hat

  po <- Y1.hat -
    Y0.hat +
    ((Y_test - Y.hat) * (W_test - W.hat)) / (W.hat * (1 - W.hat))

  list(fold = fold, po = po, Y0.hat = Y0.hat, Y1.hat = Y1.hat, W.hat = W.hat)
}

fit_tau_glmnet <- function(X, fold_list, fold_indices, po, hyper_list) {
  custom_net <- create.Learner("SL.glmnet", params = hyper_list)

  tau_list <- future_lapply(
    fold_list,
    function(fold) {
      train_filter <- fold_indices != fold
      test_filter <- !train_filter

      po_model <- SuperLearner(
        po[train_filter],
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
  fold_list
) {
  nuisances <- future_lapply(
    fold_list,
    crossfit_single_SL,
    Y,
    W,
    X,
    hyper_list,
    fold_indices,
    future.seed = TRUE
  )

  matrices <- scatter_folds(
    nuisances, fold_indices, c("po", "Y0.hat", "Y1.hat", "W.hat")
  )

  tau <- fit_tau_SL(X, fold_list, fold_indices, matrices$po, hyper_list)

  c(list(tau = tau), matrices)
}

crossfit_single_SL <- function(fold, Y, W, X, hyper_list, fold_indices) {
  train_filter <- fold_indices != fold
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
  Y0.hat <- as.numeric(predict(Y_model, mutate(test_data, W = 0L))$pred)
  Y1.hat <- as.numeric(predict(Y_model, mutate(test_data, W = 1L))$pred)
  W.hat <- trim_ps(as.numeric(predict(W_model, X[test_filter, ])$pred))

  W_test <- W[test_filter]
  Y_test <- Y[test_filter]
  Y.hat <- W_test * Y1.hat + (1 - W_test) * Y0.hat

  po <- Y1.hat -
    Y0.hat +
    ((Y_test - Y.hat) * (W_test - W.hat)) / (W.hat * (1 - W.hat))

  list(fold = fold, po = po, Y0.hat = Y0.hat, Y1.hat = Y1.hat, W.hat = W.hat)
}

fit_tau_SL <- function(X, fold_list, fold_indices, po, hyper_list) {
  tau_list <- future_lapply(
    fold_list,
    function(fold) {
      train_filter <- fold_indices != fold
      test_filter <- !train_filter

      hyper_po <- hyper_list
      hyper_po$family <- gaussian()
      hyper_po$X <- X[train_filter, ]
      hyper_po$Y <- po[train_filter]

      po_model <- do.call(SuperLearner, hyper_po)
      predict(po_model, X[test_filter, ])$pred
    },
    future.seed = TRUE
  )

  collect_tau(tau_list, fold_list, fold_indices)
}

#' The 9 candidate configurations, keyed by name
#'
#' Extracted from run_all_candidate_models() so the single-split entry point
#' below cannot drift from the crossfit one. The two differ in HOW each
#' candidate is fit, and must not differ in WHAT is being fit - the 80:20 arm
#' is a comparison of fitting strategies, so a hyperparameter that disagreed
#' between the two would confound exactly the thing it is measuring.
#'
#' @param p ncol(X)
candidate_hyperparams <- function(p) {
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
  list(
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
}

#' Which family a candidate name belongs to
#'
#' Shared by both entry points, so a new candidate needs adding in one place.
candidate_family <- function(nm) {
  if (grepl("^rf", nm)) "rf" else if (grepl("^net", nm)) "net" else "SL"
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
run_all_candidate_models <- function(Y, W, X, fold_indices, fold_list) {
  hyperparams <- candidate_hyperparams(ncol(X))

  fitter_for <- function(nm) {
    switch(candidate_family(nm), rf = fit_rf, net = fit_glmnet, SL = fit_SL)
  }

  out <- lapply(names(hyperparams), function(nm) {
    fitter_for(nm)(Y, W, X, hyper_list = hyperparams[[nm]], fold_indices, fold_list)
  })

  setNames(out, names(hyperparams))
}

# ---- single-split fitting, for the 80:20 arm --------------------------------
# The one place in this study where the candidates are fit differently. The
# 80:20 arm exists to decouple the evaluation nuisance from the candidate
# entirely: the candidate sees only the 80%, the nuisance sees only the 20%,
# so neither can borrow information from the other's rows. That is only true
# if the candidate genuinely never touches the 20% - which means refitting,
# not reusing me_analysis.R's crossfit tau (whose training data covers all n).
#
# WHAT IS AND IS NOT DIFFERENT. The estimator family is unchanged: each
# candidate is still a DR-learner, still built from leave-one-fold-out
# nuisances feeding a stage-2 regression, still with the hyperparameters
# candidate_hyperparams() defines. What changes is that the crossfit happens
# WITHIN the 80% - producing pseudo-outcomes for the training rows only - and
# stage 2 is then fit once on all of them and predicted on the untouched 20%,
# rather than fold-wise back onto the training rows themselves. A fold-wise
# stage 2 has nowhere to send its predictions here; the 20% is not one of its
# folds.

#' One candidate, fit on `train` and predicted on `test`
#'
#' @param crossfit_fn the family's per-fold nuisance function
#'   (crossfit_single_RF / _glmnet / _SL)
#' @param stage2_fn the family's whole-training-set stage-2 fitter
#' @return list(tau = predictions on `test`, plus the training-row nuisances,
#'   named *_train so nothing downstream can mistake them for length-n vectors)
fit_split_generic <- function(
  Y, W, X, hyper_list, train, test, n_folds,
  crossfit_fn, stage2_fn
) {
  Y_tr <- Y[train]
  W_tr <- W[train]
  X_tr <- X[train, , drop = FALSE]

  folds <- split_folds(Y_tr, k = n_folds)

  nuisances <- future_lapply(
    folds$fold_list,
    crossfit_fn,
    Y_tr,
    W_tr,
    X_tr,
    hyper_list,
    folds$fold_indices,
    future.seed = TRUE
  )

  matrices <- scatter_folds(
    nuisances, folds$fold_indices, c("po", "Y0.hat", "Y1.hat", "W.hat")
  )

  tau <- stage2_fn(X_tr, matrices$po, X[test, , drop = FALSE], hyper_list)

  list(
    tau = as.numeric(tau),
    po_train = matrices$po,
    Y0.hat_train = matrices$Y0.hat,
    Y1.hat_train = matrices$Y1.hat,
    W.hat_train = matrices$W.hat,
    fold_indices_train = folds$fold_indices
  )
}

stage2_split_rf <- function(X_train, po, X_test, hyper_list) {
  hyper <- hyper_list
  hyper$write.forest <- TRUE
  hyper$data <- data.frame(po = po, X_train)
  hyper$formula <- po ~ .
  predict(do.call(ranger, hyper), X_test)$predictions
}

stage2_split_glmnet <- function(X_train, po, X_test, hyper_list) {
  custom_net <- create.Learner("SL.glmnet", params = hyper_list)
  po_model <- SuperLearner(
    po,
    X_train,
    family = gaussian(),
    SL.library = custom_net$names,
    method = "method.CC_LS"
  )
  predict(po_model, X_test)$pred
}

stage2_split_SL <- function(X_train, po, X_test, hyper_list) {
  hyper <- hyper_list
  hyper$family <- gaussian()
  hyper$X <- X_train
  hyper$Y <- po
  po_model <- do.call(SuperLearner, hyper)
  predict(po_model, X_test)$pred
}

#' Fit all 9 candidates on `train`, predict on `test`
#'
#' The single-split counterpart of run_all_candidate_models(), sharing its
#' hyperparameters and its family dispatch.
#'
#' @param train,test row indices, as from single_split()
#' @return named list, one element per candidate, `tau` of length(test)
run_all_candidate_models_split <- function(
  Y, W, X, train, test, n_folds = 10L
) {
  hyperparams <- candidate_hyperparams(ncol(X))

  parts_for <- function(nm) {
    switch(
      candidate_family(nm),
      rf  = list(crossfit_single_RF, stage2_split_rf),
      net = list(crossfit_single_glmnet, stage2_split_glmnet),
      SL  = list(crossfit_single_SL, stage2_split_SL)
    )
  }

  out <- lapply(names(hyperparams), function(nm) {
    parts <- parts_for(nm)
    hyper <- hyperparams[[nm]]
    # fit_rf() sets this before crossfitting; the RF family needs it here too
    # or ranger returns a forest that cannot predict.
    if (candidate_family(nm) == "rf") hyper$write.forest <- TRUE

    fit_split_generic(
      Y, W, X, hyper, train, test, n_folds,
      crossfit_fn = parts[[1]], stage2_fn = parts[[2]]
    )
  })

  setNames(out, names(hyperparams))
}
