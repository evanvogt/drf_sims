##########
# title: crossfitting comparison - CATE model variants
##########
# Compares double crossfitting (the procedure used throughout this study) against
# standard crossfitting with the final model fit either on the whole dataset or
# through a second crossfitting pass, and - for forests - against out-of-bag
# predictions with no sample splitting at all.
#
# Stage 1 = nuisance / pseudo-outcome construction, stage 2 = final CATE regression.
# Nuisances are computed once by run_all_crossfit_variants and shared across the
# arms that use them, so per-arm timings decompose into nuisance + stage 2.

require(grf)
require(SuperLearner)
require(future)
require(furrr)
require(dplyr)
require(here)

# reused rather than forked: collate_predictions (R/utils.R), the continuous DGP
# (cts_dgms.R), and pretest_superlearner plus the reference DR implementations
# that cf_testing.R's regression check compares against (R/cate_models.R - these
# used to live in continuous/cts_models.R, which is now a thin profile shim)
source(here("R", "utils.R"))
source(here("continuous", "cts_dgms.R"))
source(here("R", "cate_models.R"))

# ---- data -------------------------------------------------------------------

#' Draw a training sample and an independent test sample from the same DGP
#'
#' bW is derived from n inside generate_continuous_scenario_data (cts_dgms.R:67),
#' so the test sample is built by stacking n_test/n draws at the same n rather
#' than by asking for a bigger sample - that keeps the true CATE surface identical.
#'
#' @param scenario scenario index passed to generate_continuous_scenario_data
#' @param n training sample size
#' @param n_test test sample size; must be a multiple of n
generate_cf_replicate <- function(scenario, n, n_test = 2000) {
  stopifnot(n_test %% n == 0)

  gen <- generate_continuous_scenario_data(scenario, n)

  reps <- lapply(seq_len(n_test / n), function(i) generate_continuous_scenario_data(scenario, n))
  test_data <- do.call(rbind, lapply(reps, `[[`, "dataset"))
  test_truth <- do.call(rbind, lapply(reps, `[[`, "truth"))

  stopifnot(all(vapply(reps, `[[`, numeric(1), "bW") == gen$bW))

  list(data = gen$dataset,
       truth_tau = gen$truth$tau,
       X_test = as.matrix(test_data[, -c(1:2)]),
       truth_test_tau = test_truth$tau,
       bW = gen$bW)
}

# ---- helpers ----------------------------------------------------------------

# AIPW / DR pseudo outcome
dr_pseudo <- function(Y, W, Y1.hat, Y0.hat, W.hat) {
  Y.hat <- W * Y1.hat + (1 - W) * Y0.hat
  cate <- Y1.hat - Y0.hat
  cate + ((Y - Y.hat) * (W - W.hat)) / (W.hat * (1 - W.hat))
}

# trim propensities away from 0/1. cts_models.R only does this for SuperLearner;
# here every arm is trimmed identically so that no arm blows up on a technicality.
# W is randomised 0.5 in the DGP so this binds only for the in-sample nuisances.
trim_ps <- function(p, lo = 0.05, hi = 0.95) pmin(pmax(p, lo), hi)

# evaluate an expression, returning its value alongside elapsed seconds
timed <- function(expr) {
  t <- system.time(val <- expr)
  list(value = val, time = unname(t["elapsed"]))
}

# reassemble leave-one-fold-out predictions into full length vectors
scatter_folds <- function(reslist, fold_indices, targets) {
  out <- lapply(targets, function(nm) {
    v <- rep(NA_real_, length(fold_indices))
    for (res in reslist) v[fold_indices == res$fold] <- res[[nm]]
    v
  })
  names(out) <- targets
  out
}

# maintainer-endorsed (unsupported) workaround for the fact that grf has no public
# API for an OOB prediction at a counterfactual covariate row: predict(forest)
# without newdata re-reads object$X.orig fresh each call and restricts each row to
# its out-of-bag trees, so pointing X.orig at a perturbed matrix and clearing the
# cached predictions borrows grf's own OOB routine at that perturbed point.
# "Hacky," per grf-labs/grf#307 - may break on a future grf version. R's
# copy-on-modify semantics mean this cannot leak into the caller's forest object;
# only the local copy inside this function is touched.
oob_predict_counterfactual <- function(forest, X_counterfactual) {
  forest$X.orig <- X_counterfactual
  forest$predictions <- NULL
  forest$debiased.error <- NULL
  predict(forest)$predictions
}

# fully-manual alternative to oob_predict_counterfactual, touching only grf's
# documented get_tree()/get_leaf_node() API rather than the X.orig monkey-patch - a
# check on whether the shortcut actually reproduces grf's own OOB routine. For each
# tree, row i is OOB if excluded from that tree's full in-bag set; the leaf a
# counterfactual point falls into is found with get_leaf_node, and predicted as the
# mean training Y over that leaf's samples. OOB status MUST be read from
# tree$drawn_samples (the tree's complete in-bag draw) rather than a leaf's own
# $samples, which for an honest forest is only the honesty-estimation subsample
# (J2) and would misclassify honesty-split-only rows as OOB (grf-labs/grf#258).
# Slower than oob_predict_counterfactual by construction (R-level loop over
# num.trees); this is a correctness check, not meant for routine use.
oob_predict_counterfactual_manual <- function(forest, X_counterfactual, Y) {
  n <- nrow(forest$X.orig)
  num_trees <- forest[["_num_trees"]]

  pred_sum <- numeric(n)
  pred_n <- integer(n)

  for (k in seq_len(num_trees)) {
    tree <- get_tree(forest, k)
    in_bag <- unique(tree$drawn_samples)
    oob_rows <- setdiff(seq_len(n), in_bag)
    if (length(oob_rows) == 0) next

    leaf_ids <- get_leaf_node(tree, X_counterfactual[oob_rows, , drop = FALSE])
    for (j in seq_along(oob_rows)) {
      leaf_samples <- tree$nodes[[leaf_ids[j]]]$samples
      pred_sum[oob_rows[j]] <- pred_sum[oob_rows[j]] + mean(Y[leaf_samples])
      pred_n[oob_rows[j]] <- pred_n[oob_rows[j]] + 1
    }
  }

  ifelse(pred_n > 0, pred_sum / pred_n, NA_real_)
}

# ---- stage 1: nuisance estimation (random forest) ---------------------------

# double crossfitting over fold pairs - the status quo (cts_models.R:62)
nuisance_double_rf <- function(X, Y, W, fold_indices, fold_pairs, num.threads = NULL) {

  cross_fits <- future_map(seq_along(fold_pairs), function(i) {
    fold_pair <- fold_pairs[[i]]
    in_train <- !(fold_indices %in% fold_pair)
    in_test <- !in_train

    Y.hat.model <- regression_forest(cbind(W = W[in_train], X[in_train, ]), Y[in_train],
                                     num.threads = num.threads)
    Y.hat.cf.model <- regression_forest(X[in_train, ], Y[in_train], num.threads = num.threads)
    W.hat.model <- regression_forest(X[in_train, ], W[in_train], num.threads = num.threads)

    X_test <- X[in_test, ]
    Y0.hat <- predict(Y.hat.model, newdata = cbind(W = 0, X_test))$predictions
    Y1.hat <- predict(Y.hat.model, newdata = cbind(W = 1, X_test))$predictions
    Y.hat.cf <- predict(Y.hat.cf.model, newdata = X_test)$predictions
    W.hat <- trim_ps(predict(W.hat.model, newdata = X_test)$predictions)

    list(po = dr_pseudo(Y[in_test], W[in_test], Y1.hat, Y0.hat, W.hat),
         Y0.hat = Y0.hat, Y.hat.cf = Y.hat.cf, W.hat = W.hat)
  }, .options = furrr_options(seed = TRUE))

  fold_list <- unique(fold_indices)
  po_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "po")
  Y.hat.cf_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "Y.hat.cf")
  Y0.hat_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "Y0.hat")
  W.hat_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "W.hat")

  list(po = po_matrix,
       Y.hat.cf_matrix = Y.hat.cf_matrix,
       W.hat_matrix = W.hat_matrix,
       Y0.hat = rowMeans(Y0.hat_matrix, na.rm = TRUE),
       Y.hat.cf = rowMeans(Y.hat.cf_matrix, na.rm = TRUE),
       W.hat = rowMeans(W.hat_matrix, na.rm = TRUE))
}

# ordinary leave-one-fold-out crossfitting, S-learner outcome model
nuisance_single_rf <- function(X, Y, W, fold_indices, num.threads = NULL) {

  cross_fits <- future_map(unique(fold_indices), function(fold) {
    in_train <- fold_indices != fold
    in_test <- !in_train

    Y.hat.model <- regression_forest(cbind(W = W[in_train], X[in_train, ]), Y[in_train],
                                     num.threads = num.threads)
    Y.hat.cf.model <- regression_forest(X[in_train, ], Y[in_train], num.threads = num.threads)
    W.hat.model <- regression_forest(X[in_train, ], W[in_train], num.threads = num.threads)

    X_test <- X[in_test, ]
    Y0.hat <- predict(Y.hat.model, newdata = cbind(W = 0, X_test))$predictions
    Y1.hat <- predict(Y.hat.model, newdata = cbind(W = 1, X_test))$predictions
    Y.hat.cf <- predict(Y.hat.cf.model, newdata = X_test)$predictions
    W.hat <- trim_ps(predict(W.hat.model, newdata = X_test)$predictions)

    list(fold = fold, po = dr_pseudo(Y[in_test], W[in_test], Y1.hat, Y0.hat, W.hat),
         Y0.hat = Y0.hat, Y.hat.cf = Y.hat.cf, W.hat = W.hat)
  }, .options = furrr_options(seed = TRUE))

  scatter_folds(cross_fits, fold_indices, c("po", "Y0.hat", "Y.hat.cf", "W.hat"))
}

# leave-one-fold-out crossfitting with a T-learner outcome model. exists as the
# control for oob_oob: it differs from nuisance_single_rf only in learner
# structure, so the S-vs-T effect and the OOB-vs-crossfit effect are separable.
nuisance_single_rf_t <- function(X, Y, W, fold_indices, num.threads = NULL) {

  cross_fits <- future_map(unique(fold_indices), function(fold) {
    in_train <- fold_indices != fold
    in_test <- !in_train

    X_train <- X[in_train, , drop = FALSE]
    Y_train <- Y[in_train]
    W_train <- W[in_train]

    f1 <- regression_forest(X_train[W_train == 1, , drop = FALSE], Y_train[W_train == 1],
                            num.threads = num.threads)
    f0 <- regression_forest(X_train[W_train == 0, , drop = FALSE], Y_train[W_train == 0],
                            num.threads = num.threads)
    W.hat.model <- regression_forest(X_train, W_train, num.threads = num.threads)

    X_test <- X[in_test, , drop = FALSE]
    Y1.hat <- predict(f1, newdata = X_test)$predictions
    Y0.hat <- predict(f0, newdata = X_test)$predictions
    W.hat <- trim_ps(predict(W.hat.model, newdata = X_test)$predictions)

    list(fold = fold, po = dr_pseudo(Y[in_test], W[in_test], Y1.hat, Y0.hat, W.hat),
         Y0.hat = Y0.hat, W.hat = W.hat)
  }, .options = furrr_options(seed = TRUE))

  scatter_folds(cross_fits, fold_indices, c("po", "Y0.hat", "W.hat"))
}

# no sample splitting: forests fit on the whole sample, predictions taken out-of-bag.
# a T-learner is required here because grf can only return OOB predictions at each
# unit's observed covariate row - an S-learner on cbind(W, X) has no OOB counterfactual.
# with separate arm forests, a treated unit's Y1.hat is OOB and its Y0.hat comes from
# a forest that never saw it, so both arms are honest.
nuisance_oob_rf <- function(X, Y, W, num.threads = NULL) {
  n_obs <- nrow(X)
  treated <- W == 1

  f1 <- regression_forest(X[treated, , drop = FALSE], Y[treated], num.threads = num.threads)
  f0 <- regression_forest(X[!treated, , drop = FALSE], Y[!treated], num.threads = num.threads)

  Y1.hat <- numeric(n_obs)
  Y1.hat[treated] <- predict(f1)$predictions
  Y1.hat[!treated] <- predict(f1, newdata = X[!treated, , drop = FALSE])$predictions

  Y0.hat <- numeric(n_obs)
  Y0.hat[!treated] <- predict(f0)$predictions
  Y0.hat[treated] <- predict(f0, newdata = X[treated, , drop = FALSE])$predictions

  W.hat <- trim_ps(predict(regression_forest(X, W, num.threads = num.threads))$predictions)
  Y.hat.cf <- predict(regression_forest(X, Y, num.threads = num.threads))$predictions

  list(po = dr_pseudo(Y, W, Y1.hat, Y0.hat, W.hat),
       Y0.hat = Y0.hat, Y.hat.cf = Y.hat.cf, W.hat = W.hat)
}

# S-learner counterpart to nuisance_oob_rf, using oob_predict_counterfactual instead
# of a T-learner - matches the S-learner structure used by every other nuisance in
# this file (nuisance_double_rf, nuisance_single_rf). Kept alongside nuisance_oob_rf
# rather than replacing it, so the T-learner and the workaround can be compared.
nuisance_oob_rf_s <- function(X, Y, W, num.threads = NULL) {
  forest <- regression_forest(cbind(W = W, X), Y, num.threads = num.threads)

  Y0.hat <- oob_predict_counterfactual(forest, cbind(W = 0, X))
  Y1.hat <- oob_predict_counterfactual(forest, cbind(W = 1, X))

  W.hat <- trim_ps(predict(regression_forest(X, W, num.threads = num.threads))$predictions)
  Y.hat.cf <- predict(regression_forest(X, Y, num.threads = num.threads))$predictions

  list(po = dr_pseudo(Y, W, Y1.hat, Y0.hat, W.hat),
       Y0.hat = Y0.hat, Y.hat.cf = Y.hat.cf, W.hat = W.hat)
}

# manual-API counterpart to nuisance_oob_rf_s - identical S-learner forest, but the
# counterfactual OOB predictions come from oob_predict_counterfactual_manual instead
# of the X.orig shortcut. Y0.hat/Y1.hat can be NA for rows that are in-bag in every
# tree (astronomically unlikely at grf's default num.trees, but dr_pseudo will
# propagate a resulting NA rather than silently drop it, so a run that hits this
# will fail loudly).
nuisance_oob_rf_manual <- function(X, Y, W, num.threads = NULL) {
  forest <- regression_forest(cbind(W = W, X), Y, num.threads = num.threads)

  Y0.hat <- oob_predict_counterfactual_manual(forest, cbind(W = 0, X), Y)
  Y1.hat <- oob_predict_counterfactual_manual(forest, cbind(W = 1, X), Y)

  W.hat <- trim_ps(predict(regression_forest(X, W, num.threads = num.threads))$predictions)
  Y.hat.cf <- predict(regression_forest(X, Y, num.threads = num.threads))$predictions

  list(po = dr_pseudo(Y, W, Y1.hat, Y0.hat, W.hat),
       Y0.hat = Y0.hat, Y.hat.cf = Y.hat.cf, W.hat = W.hat)
}

# ---- stage 1: nuisance estimation (SuperLearner) ----------------------------

# one train/test split's worth of SuperLearner nuisances. X must be a data.frame.
# failsafes and propensity trimming carried over from cts_models.R:300-314.
sl_nuisance_fit <- function(X, Y, W, in_train, in_test, sl_lib) {

  X_train <- X[in_train, , drop = FALSE]
  X_W_train <- cbind(W = W[in_train], X_train)

  Y_lib <- pretest_superlearner(Y[in_train], X_W_train, sl_lib, gaussian())
  Y.hat.model <- SuperLearner(Y = Y[in_train], X = X_W_train, SL.library = Y_lib)

  W_lib <- pretest_superlearner(W[in_train], X_train, sl_lib, binomial())
  W.hat.model <- SuperLearner(W[in_train], X_train, family = binomial(),
                              SL.library = W_lib, method = "method.NNloglik")

  X_test <- X[in_test, , drop = FALSE]
  Y0.hat <- as.numeric(predict(Y.hat.model, newdata = cbind(W = 0, X_test))$pred)
  Y1.hat <- as.numeric(predict(Y.hat.model, newdata = cbind(W = 1, X_test))$pred)
  W.hat <- as.numeric(predict(W.hat.model, newdata = X_test)$pred)

  if (all(Y0.hat == 0) && all(Y1.hat == 0)) {
    warning("SuperLearner failed for Y.hat. Using mean(Y).")
    Y0.hat <- rep(mean(Y[in_train][W[in_train] == 0], na.rm = TRUE), sum(in_test))
    Y1.hat <- rep(mean(Y[in_train][W[in_train] == 1], na.rm = TRUE), sum(in_test))
  }
  if (all(W.hat == 0)) {
    warning("SuperLearner failed for W.hat. Using mean(W).")
    W.hat <- rep(mean(W[in_train], na.rm = TRUE), sum(in_test))
  }
  W.hat <- trim_ps(W.hat)

  list(po = dr_pseudo(Y[in_test], W[in_test], Y1.hat, Y0.hat, W.hat),
       Y0.hat = Y0.hat, W.hat = W.hat)
}

nuisance_double_sl <- function(X, Y, W, fold_indices, fold_pairs, sl_lib) {

  cross_fits <- future_map(seq_along(fold_pairs), function(i) {
    fold_pair <- fold_pairs[[i]]
    in_train <- !(fold_indices %in% fold_pair)
    sl_nuisance_fit(X, Y, W, in_train, !in_train, sl_lib)
  }, .options = furrr_options(seed = TRUE))

  fold_list <- unique(fold_indices)
  po_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "po")
  Y0.hat_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "Y0.hat")
  W.hat_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "W.hat")

  list(po = po_matrix,
       Y0.hat = rowMeans(Y0.hat_matrix, na.rm = TRUE),
       W.hat = rowMeans(W.hat_matrix, na.rm = TRUE))
}

nuisance_single_sl <- function(X, Y, W, fold_indices, sl_lib) {

  cross_fits <- future_map(unique(fold_indices), function(fold) {
    in_train <- fold_indices != fold
    c(list(fold = fold), sl_nuisance_fit(X, Y, W, in_train, !in_train, sl_lib))
  }, .options = furrr_options(seed = TRUE))

  scatter_folds(cross_fits, fold_indices, c("po", "Y0.hat", "W.hat"))
}

# ---- stage 2: final CATE regression -----------------------------------------

# crossfit stage 2: fit on the complement of each fold, predict the held-out fold.
# po is either the n x V matrix from double crossfitting (column k is untouched by
# fold k) or a plain n-vector. the test-set prediction averages the V fold models.
stage2_crossfit_rf <- function(X, po, X_test, fold_indices, num.threads = NULL) {
  po_is_matrix <- is.matrix(po)
  # a po matrix is indexed by the fold it is valid for, so it is only meaningful
  # against the split it was built from - never against the fresh split
  stopifnot(!po_is_matrix || ncol(po) == length(unique(fold_indices)))

  fits <- future_map(unique(fold_indices), function(fold) {
    in_train <- fold_indices != fold
    in_fold <- !in_train
    y_train <- if (po_is_matrix) po[in_train, fold] else po[in_train]

    forest <- regression_forest(X[in_train, , drop = FALSE], y_train, num.threads = num.threads)
    list(fold = fold,
         tau = predict(forest, newdata = X[in_fold, , drop = FALSE])$predictions,
         tau_test = predict(forest, newdata = X_test)$predictions)
  }, .options = furrr_options(seed = TRUE))

  tau <- rep(NA_real_, nrow(X))
  for (fit in fits) tau[fold_indices == fit$fold] <- fit$tau

  # tau_test_folds is kept in memory only, for the single-model test score in
  # arm() - it is never saved, being n_test x V per crossfit arm
  tau_test_folds <- sapply(fits, `[[`, "tau_test")
  list(tau = tau, tau_test = rowMeans(tau_test_folds), tau_test_folds = tau_test_folds)
}

# whole-sample stage 2: one forest, its OOB predictions and its test-set predictions.
stage2_whole_rf <- function(X, po, X_test, num.threads = NULL) {
  forest <- regression_forest(X, po, num.threads = num.threads)
  list(tau_oob = predict(forest)$predictions,
       tau_test = predict(forest, newdata = X_test)$predictions)
}

stage2_crossfit_sl <- function(X, po, X_test, fold_indices, sl_lib) {
  po_is_matrix <- is.matrix(po)
  stopifnot(!po_is_matrix || ncol(po) == length(unique(fold_indices)))

  fits <- future_map(unique(fold_indices), function(fold) {
    in_train <- fold_indices != fold
    in_fold <- !in_train
    y_train <- if (po_is_matrix) po[in_train, fold] else po[in_train]
    X_train <- X[in_train, , drop = FALSE]

    # cts_models.R:387 pretests into po_lib but then passes the untested sl_lib
    # in the matrix branch; po_lib is used in both branches here
    po_lib <- pretest_superlearner(y_train, X_train, sl_lib, gaussian())
    po_model <- SuperLearner(y_train, X_train, family = gaussian(), SL.library = po_lib)

    list(fold = fold,
         tau = as.numeric(predict(po_model, newdata = X[in_fold, , drop = FALSE])$pred),
         tau_test = as.numeric(predict(po_model, newdata = X_test)$pred))
  }, .options = furrr_options(seed = TRUE))

  tau <- rep(NA_real_, nrow(X))
  for (fit in fits) tau[fold_indices == fit$fold] <- fit$tau

  tau_test_folds <- sapply(fits, `[[`, "tau_test")
  list(tau = tau, tau_test = rowMeans(tau_test_folds), tau_test_folds = tau_test_folds)
}

# ---- causal forest ----------------------------------------------------------

# fold-wise causal forest. Y.hat / W.hat are either n x V matrices indexed by fold
# (double crossfitting) or plain n-vectors (ordinary crossfitting).
cf_foldwise <- function(X, Y, W, X_test, Y.hat, W.hat, fold_indices, num.threads = NULL) {
  hat_is_matrix <- is.matrix(Y.hat)

  fits <- future_map(unique(fold_indices), function(fold) {
    in_train <- fold_indices != fold
    in_fold <- !in_train
    y_hat <- if (hat_is_matrix) Y.hat[in_train, fold] else Y.hat[in_train]
    w_hat <- if (hat_is_matrix) W.hat[in_train, fold] else W.hat[in_train]

    forest <- causal_forest(X[in_train, , drop = FALSE], Y[in_train], W[in_train],
                            y_hat, w_hat, num.threads = num.threads)
    list(fold = fold,
         tau = predict(forest, newdata = X[in_fold, , drop = FALSE])$predictions,
         tau_test = predict(forest, newdata = X_test)$predictions)
  }, .options = furrr_options(seed = TRUE))

  tau <- rep(NA_real_, nrow(X))
  for (fit in fits) tau[fold_indices == fit$fold] <- fit$tau

  tau_test_folds <- sapply(fits, `[[`, "tau_test")
  list(tau = tau, tau_test = rowMeans(tau_test_folds), tau_test_folds = tau_test_folds)
}

# whole-sample causal forest: its OOB predictions and its test-set predictions.
# Y.hat / W.hat NULL falls back to grf's own internally cross-fit OOB nuisances.
cf_whole <- function(X, Y, W, X_test, Y.hat = NULL, W.hat = NULL, num.threads = NULL) {
  forest <- causal_forest(X, Y, W, Y.hat = Y.hat, W.hat = W.hat, num.threads = num.threads)
  list(tau_oob = predict(forest)$predictions,
       tau_test = predict(forest, newdata = X_test)$predictions)
}

# ---- orchestrator -----------------------------------------------------------

# assemble one arm's record.
#
# mse_test_single exists because the crossfit and whole-sample arms do not predict
# on the test set in the same way. A crossfit arm ends up with V fitted models and
# averages their test predictions - that is the estimator you would actually deploy,
# so tau_test keeps it - but the averaging is a variance-reducing ensemble on top of
# the honesty effect being studied. mse_test_single scores each fold model on the
# test set separately and averages the V scores, which is the like-for-like reading
# against a whole-sample arm's single model. For whole-sample arms the two coincide.
arm <- function(family, variant, tau, tau_test, time_nuisance, time_stage2,
                tau_test_folds = NULL, truth_test = NULL) {

  mse_test_single <- if (is.null(truth_test)) {
    NA_real_
  } else if (is.null(tau_test_folds)) {
    mean((tau_test - truth_test)^2)
  } else {
    mean(apply(tau_test_folds, 2, function(p) mean((p - truth_test)^2)))
  }

  list(family = family, variant = variant, tau = tau, tau_test = tau_test,
       time_nuisance = time_nuisance, time_stage2 = time_stage2,
       mse_test_single = mse_test_single)
}

#' Run every crossfitting variant on one simulated dataset
#'
#' @param data data.frame with columns Y, W then covariates (the DGM layout)
#' @param X_test covariate matrix for the independent test sample
#' @param n_folds number of folds V
#' @param sl_lib SuperLearner library; NULL skips the SuperLearner family
#' @param num.threads grf thread count; NULL is grf's default (all cores)
#' @param truth_test true CATEs on the test sample. Used for scoring only - it is
#'   never seen by any model - and only to compute the per-arm mse_test_single.
#'   NULL leaves that field NA.
#' @return list with $arms (named list of arm records), $fold_indices, $fold_indices_b
run_all_crossfit_variants <- function(data, X_test, n_folds = 10, sl_lib = NULL,
                                      num.threads = NULL, truth_test = NULL) {

  X <- as.matrix(data[, -c(1:2)])
  Y <- data$Y
  W <- data$W
  n_obs <- nrow(X)

  # one fold assignment shared by every variant, so differences are attributable
  # to the procedure and not to the fold draw. rows are i.i.d. so the deterministic
  # assignment used throughout this study is fine.
  fold_indices <- sort(seq(n_obs) %% n_folds) + 1
  fold_list <- unique(fold_indices)
  fold_pairs <- utils::combn(fold_list, 2, simplify = FALSE)

  # a fresh, independent split for the scf_scf_new variant
  fold_indices_b <- sample(fold_indices)

  arms <- list()

  # -- nuisances (random forest) ----------------------------------------------
  cat("Nuisances: RF double crossfit...\n")
  nz_double <- timed(nuisance_double_rf(X, Y, W, fold_indices, fold_pairs, num.threads))
  cat("Nuisances: RF single crossfit...\n")
  nz_single <- timed(nuisance_single_rf(X, Y, W, fold_indices, num.threads))
  cat("Nuisances: RF single crossfit (T-learner)...\n")
  nz_single_t <- timed(nuisance_single_rf_t(X, Y, W, fold_indices, num.threads))
  cat("Nuisances: RF out-of-bag...\n")
  nz_oob <- timed(nuisance_oob_rf(X, Y, W, num.threads))
  cat("Nuisances: RF out-of-bag (S-learner, X.orig workaround)...\n")
  nz_oob_s <- timed(nuisance_oob_rf_s(X, Y, W, num.threads))
  cat("Nuisances: RF out-of-bag (S-learner, manual tree-loop)...\n")
  nz_oob_manual <- timed(nuisance_oob_rf_manual(X, Y, W, num.threads))

  # -- DR learner, random forest ---------------------------------------------
  cat("DR-RF variants...\n")

  s <- timed(stage2_crossfit_rf(X, nz_double$value$po, X_test, fold_indices, num.threads))
  arms$dcf <- arm("dr_rf", "dcf", s$value$tau, s$value$tau_test, nz_double$time, s$time,
                  s$value$tau_test_folds, truth_test)

  s <- timed(stage2_crossfit_rf(X, nz_single$value$po, X_test, fold_indices, num.threads))
  arms$scf_scf <- arm("dr_rf", "scf_scf", s$value$tau, s$value$tau_test, nz_single$time, s$time,
                      s$value$tau_test_folds, truth_test)

  s <- timed(stage2_crossfit_rf(X, nz_single$value$po, X_test, fold_indices_b, num.threads))
  arms$scf_scf_new <- arm("dr_rf", "scf_scf_new", s$value$tau, s$value$tau_test,
                          nz_single$time, s$time, s$value$tau_test_folds, truth_test)

  s <- timed(stage2_whole_rf(X, nz_single$value$po, X_test, num.threads))
  arms$scf_oob <- arm("dr_rf", "scf_oob", s$value$tau_oob, s$value$tau_test,
                      nz_single$time, s$time, truth_test = truth_test)

  s <- timed(stage2_whole_rf(X, nz_single_t$value$po, X_test, num.threads))
  arms$scf_oob_t <- arm("dr_rf", "scf_oob_t", s$value$tau_oob, s$value$tau_test,
                        nz_single_t$time, s$time, truth_test = truth_test)

  s <- timed(stage2_whole_rf(X, nz_oob$value$po, X_test, num.threads))
  arms$oob_oob <- arm("dr_rf", "oob_oob", s$value$tau_oob, s$value$tau_test,
                      nz_oob$time, s$time, truth_test = truth_test)

  s <- timed(stage2_whole_rf(X, nz_oob_s$value$po, X_test, num.threads))
  arms$oob_oob_s <- arm("dr_rf", "oob_oob_s", s$value$tau_oob, s$value$tau_test,
                        nz_oob_s$time, s$time, truth_test = truth_test)

  s <- timed(stage2_whole_rf(X, nz_oob_manual$value$po, X_test, num.threads))
  arms$oob_oob_manual <- arm("dr_rf", "oob_oob_manual", s$value$tau_oob, s$value$tau_test,
                             nz_oob_manual$time, s$time, truth_test = truth_test)

  # -- causal forest ---------------------------------------------------------
  cat("Causal forest variants...\n")

  s <- timed(cf_foldwise(X, Y, W, X_test, nz_double$value$Y.hat.cf_matrix,
                         nz_double$value$W.hat_matrix, fold_indices, num.threads))
  arms$cf_dcf <- arm("causal_forest", "cf_dcf", s$value$tau, s$value$tau_test,
                     nz_double$time, s$time, s$value$tau_test_folds, truth_test)

  s <- timed(cf_foldwise(X, Y, W, X_test, nz_single$value$Y.hat.cf,
                         nz_single$value$W.hat, fold_indices, num.threads))
  arms$cf_scf <- arm("causal_forest", "cf_scf", s$value$tau, s$value$tau_test,
                     nz_single$time, s$time, s$value$tau_test_folds, truth_test)

  s <- timed(cf_whole(X, Y, W, X_test, nz_single$value$Y.hat.cf,
                      nz_single$value$W.hat, num.threads))
  arms$cf_full_oob <- arm("causal_forest", "cf_full_oob", s$value$tau_oob, s$value$tau_test,
                          nz_single$time, s$time, truth_test = truth_test)

  # grf's own internally cross-fit OOB nuisances - no separate nuisance stage
  s <- timed(cf_whole(X, Y, W, X_test, num.threads = num.threads))
  arms$cf_default <- arm("causal_forest", "cf_default", s$value$tau_oob, s$value$tau_test,
                         0, s$time, truth_test = truth_test)

  # -- DR learner, SuperLearner ----------------------------------------------
  # no OOB analogue exists for SuperLearner, so the OOB arms and the T-learner
  # control are dropped from this family
  if (!is.null(sl_lib)) {
    X_df <- as.data.frame(X)
    X_test_df <- as.data.frame(X_test)

    cat("Nuisances: SL double crossfit...\n")
    sz_double <- timed(nuisance_double_sl(X_df, Y, W, fold_indices, fold_pairs, sl_lib))
    cat("Nuisances: SL single crossfit...\n")
    sz_single <- timed(nuisance_single_sl(X_df, Y, W, fold_indices, sl_lib))

    cat("DR-SL variants...\n")

    s <- timed(stage2_crossfit_sl(X_df, sz_double$value$po, X_test_df, fold_indices, sl_lib))
    arms$sl_dcf <- arm("dr_sl", "dcf", s$value$tau, s$value$tau_test, sz_double$time, s$time,
                       s$value$tau_test_folds, truth_test)

    s <- timed(stage2_crossfit_sl(X_df, sz_single$value$po, X_test_df, fold_indices, sl_lib))
    arms$sl_scf_scf <- arm("dr_sl", "scf_scf", s$value$tau, s$value$tau_test,
                           sz_single$time, s$time, s$value$tau_test_folds, truth_test)

    s <- timed(stage2_crossfit_sl(X_df, sz_single$value$po, X_test_df, fold_indices_b, sl_lib))
    arms$sl_scf_scf_new <- arm("dr_sl", "scf_scf_new", s$value$tau, s$value$tau_test,
                               sz_single$time, s$time, s$value$tau_test_folds, truth_test)
  }

  list(arms = arms, fold_indices = fold_indices, fold_indices_b = fold_indices_b)
}
