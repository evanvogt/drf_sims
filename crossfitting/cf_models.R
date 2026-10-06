##########
# title: crossfitting comparison - CATE model variants
##########
# Compares double crossfitting (the procedure this study used to use throughout)
# against standard crossfitting with the final model fit either on the whole
# dataset or through a second crossfitting pass, and - for forests - against
# out-of-bag predictions with no sample splitting at all.
#
# Every DR arm fits its outcome model separately in each treatment arm (a
# T-learner), as the production DR-learners in R/cate_models.R do, so the arms
# differ only in how they split the sample. The per-arm fits themselves are
# production's: t_learner_rf_split / nuisance_rf for the forests, sl_split_fit /
# nuisance_sl for SuperLearner.
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
# (the "continuous" set of R/dgm_scenarios.R), and from R/cate_models.R the
# per-arm outcome models (t_learner_rf_split, nuisance_rf, sl_split_fit,
# nuisance_sl) plus, via R/sl_library.R, pretest_superlearner and sl_fit_predict
source(here("R", "utils.R"))
source(here("R", "dgm_scenarios.R"))
source(here("R", "cate_models.R"))

# ---- data -------------------------------------------------------------------

#' Draw a training sample and an independent test sample from the same DGP
#'
#' bW is derived from n inside generate_scenario_data (calibrate_bW() in
#' R/dgm_scenarios.R), so the test sample is built by stacking n_test/n draws
#' at the same n rather
#' than by asking for a bigger sample - that keeps the true CATE surface identical.
#'
#' @param scenario scenario index of the "continuous" set, passed to
#'   generate_scenario_data
#' @param n training sample size
#' @param n_test test sample size; must be a multiple of n
generate_cf_replicate <- function(scenario, n, n_test = 2000) {
  stopifnot(n_test %% n == 0)

  gen <- generate_scenario_data(scenario, n, set = "continuous")

  reps <- lapply(seq_len(n_test / n),
                 function(i) generate_scenario_data(scenario, n, set = "continuous"))
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

# trim propensities away from 0/1, in every arm identically so that no arm blows
# up on a technicality (R/cate_models.R now trims the same way, RF and SL alike).
# W is randomised 0.5 in the DGP so this rarely binds.
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

# ---- stage 1: nuisance estimation (random forest) ---------------------------

# one train/test split's RF nuisances. The outcome model is production's per-arm
# forest pair (R/cate_models.R's t_learner_rf_split, control arm first);
# Y.hat.cf is the marginal outcome forest the causal forest arms need, and W.hat
# the propensity. Shared by the double and single crossfit arms, so they differ
# only in which rows each split trains on.
rf_split_fit <- function(X, Y, W, in_train, in_test, num.threads = NULL) {
  mu <- t_learner_rf_split(X, Y, W, in_train, in_test, num.threads = num.threads)

  X_train <- X[in_train, , drop = FALSE]
  X_test <- X[in_test, , drop = FALSE]
  Y.hat.cf.model <- regression_forest(X_train, Y[in_train], num.threads = num.threads)
  W.hat.model <- regression_forest(X_train, W[in_train], num.threads = num.threads)
  Y.hat.cf <- predict(Y.hat.cf.model, newdata = X_test)$predictions
  W.hat <- trim_ps(predict(W.hat.model, newdata = X_test)$predictions)

  list(po = dr_pseudo(Y[in_test], W[in_test], mu$Y1.hat, mu$Y0.hat, W.hat),
       Y0.hat = mu$Y0.hat, Y.hat.cf = Y.hat.cf, W.hat = W.hat)
}

# double crossfitting over fold pairs - the status quo this study benchmarks
nuisance_double_rf <- function(X, Y, W, fold_indices, fold_pairs, num.threads = NULL) {

  cross_fits <- future_map(seq_along(fold_pairs), function(i) {
    in_train <- !(fold_indices %in% fold_pairs[[i]])
    rf_split_fit(X, Y, W, in_train, !in_train, num.threads)
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

# ordinary leave-one-fold-out crossfitting
nuisance_single_rf <- function(X, Y, W, fold_indices, num.threads = NULL) {

  cross_fits <- future_map(unique(fold_indices), function(fold) {
    in_train <- fold_indices != fold
    c(list(fold = fold), rf_split_fit(X, Y, W, in_train, !in_train, num.threads))
  }, .options = furrr_options(seed = TRUE))

  scatter_folds(cross_fits, fold_indices, c("po", "Y0.hat", "Y.hat.cf", "W.hat"))
}

# no sample splitting: production's own whole-sample nuisances (R/cate_models.R's
# nuisance_rf), so the oob_oob arm is dr_random_forest's stage 1 by construction.
# Per-arm forests make OOB honest without any workaround: a unit's own-arm
# prediction is OOB and its other-arm prediction comes from a forest that never
# saw it.
nuisance_oob_rf <- function(X, Y, W, num.threads = NULL) {
  nuisance_rf(X, Y, W, num.threads = num.threads)
}

# ---- stage 1: nuisance estimation (SuperLearner) ----------------------------

# X must be a data.frame. Each split's fit is production's sl_split_fit
# (R/cate_models.R) - one SuperLearner per arm, plus the propensity - so the
# dcf arm differs from scf_scf only in splitting. sl_lib is list(W = , Y = ,
# tau = ) or one character vector - see R/sl_library.R's as_sl_libs(). Dropped
# learners are not recorded here, as they never were for this arm.
nuisance_double_sl <- function(X, Y, W, fold_indices, fold_pairs, sl_lib) {

  sl_lib <- as_sl_libs(sl_lib)

  cross_fits <- future_map(seq_along(fold_pairs), function(i) {
    in_train <- !(fold_indices %in% fold_pairs[[i]])
    sl_split_fit(X, Y, W, in_train, !in_train, sl_lib)
  }, .options = furrr_options(seed = TRUE))

  fold_list <- unique(fold_indices)
  po_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "po")
  Y0.hat_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "Y0.hat")
  W.hat_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "W.hat")

  list(po = po_matrix,
       Y0.hat = rowMeans(Y0.hat_matrix, na.rm = TRUE),
       W.hat = rowMeans(W.hat_matrix, na.rm = TRUE))
}

# single crossfit: production's nuisance_sl itself, so the scf_scf arm's stage 1
# is dr_superlearner's by construction
nuisance_single_sl <- function(X, Y, W, fold_indices, sl_lib) {
  nuisance_sl(X, Y, W, fold_indices, sl_lib)
}

# ---- stage 2: final CATE regression -----------------------------------------

# crossfit stage 2: fit on the complement of each fold, predict the held-out fold.
# po is either the n x V matrix from double crossfitting (column k is untouched by
# fold k) or a plain n-vector. the test-set prediction averages the V fold models.
stage2_crossfit_rf <- function(X, po, X_test, fold_indices, num.threads = NULL) {
  po_is_matrix <- is.matrix(po)
  # a po matrix is indexed by the fold it is valid for, so it is only meaningful
  # against the split it was built from
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
#
# var_oob is grf's own OOB variance estimate (bootstrap of little bags), free
# alongside the predictions since regression_forest already defaults to
# ci.group.size = 2. It gives the OOB arms a second, closed-form interval to
# score against the half-sample bootstrap band - see R/metrics.R's
# normal_interval. Note it treats po as a known outcome, so like the bootstrap it
# carries no first-stage nuisance uncertainty. predict() does not consume R's RNG
# stream, so asking for the variance leaves every point estimate unchanged.
stage2_whole_rf <- function(X, po, X_test, num.threads = NULL) {
  forest <- regression_forest(X, po, num.threads = num.threads)
  p_oob <- predict(forest, estimate.variance = TRUE)
  list(tau_oob = p_oob$predictions,
       var_oob = p_oob$variance.estimates,
       tau_test = predict(forest, newdata = X_test)$predictions)
}

stage2_crossfit_sl <- function(X, po, X_test, fold_indices, sl_lib) {
  po_is_matrix <- is.matrix(po)
  stopifnot(!po_is_matrix || ncol(po) == length(unique(fold_indices)))
  tau_lib <- as_sl_libs(sl_lib)$tau

  fits <- future_map(unique(fold_indices), function(fold) {
    in_train <- fold_indices != fold
    in_fold <- !in_train
    y_train <- if (po_is_matrix) po[in_train, fold] else po[in_train]
    X_train <- X[in_train, , drop = FALSE]

    # cts_models.R:387 pretests into po_lib but then passes the untested sl_lib
    # in the matrix branch; po_lib is used in both branches here
    po_lib <- pretest_superlearner(y_train, X_train, tau_lib, gaussian())
    po_fit <- sl_fit_predict(y_train, X_train,
                             list(fold = X[in_fold, , drop = FALSE], test = X_test),
                             po_lib)

    list(fold = fold, tau = po_fit$pred$fold, tau_test = po_fit$pred$test)
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
#
# var_oob mirrors stage2_whole_rf's. The fitted hats are returned too, named to
# the nuisance-object convention: cf_default supplies none of its own, so this is
# the only way its half-sample bootstrap can hold nuisances fixed the way every
# other arm's does (R/bootstrap_ci.R's cf_oob_half_boot).
cf_whole <- function(X, Y, W, X_test, Y.hat = NULL, W.hat = NULL, num.threads = NULL) {
  forest <- causal_forest(X, Y, W, Y.hat = Y.hat, W.hat = W.hat, num.threads = num.threads)
  p_oob <- predict(forest, estimate.variance = TRUE)
  list(tau_oob = p_oob$predictions,
       var_oob = p_oob$variance.estimates,
       tau_test = predict(forest, newdata = X_test)$predictions,
       Y.hat.cf = forest$Y.hat, W.hat = forest$W.hat)
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
#
# variance is grf's own OOB variance estimate, and only the whole-sample/OOB arms
# have one - a crossfit arm's tau is stitched together from V forests each
# predicting its own held-out fold, which is not the quantity grf's variance
# theory covers. NULL for those arms, and downstream code keys off that: see
# crossfitting/confidence_intervals/cf_ci_metrics.R, which emits a grf_normal
# interval row exactly where this field is present.
arm <- function(family, variant, tau, tau_test, time_nuisance, time_stage2,
                tau_test_folds = NULL, truth_test = NULL, variance = NULL) {

  mse_test_single <- if (is.null(truth_test)) {
    NA_real_
  } else if (is.null(tau_test_folds)) {
    mean((tau_test - truth_test)^2)
  } else {
    mean(apply(tau_test_folds, 2, function(p) mean((p - truth_test)^2)))
  }

  list(family = family, variant = variant, tau = tau, tau_test = tau_test,
       time_nuisance = time_nuisance, time_stage2 = time_stage2,
       mse_test_single = mse_test_single, var_oob = variance)
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
#' @return list with $arms (named list of arm records), $fold_indices, and
#'   $nuisances (the raw nuisance objects, needed by the
#'   half-sample bootstraps in R/bootstrap_ci.R - see
#'   crossfitting/confidence_intervals/cf_ci_analysis.R, which calls this with
#'   sl_lib = NULL and bootstraps all 8 RF/CF arms). $nuisances is large; every
#'   caller picks the fields it saves, and cf_analysis.R deliberately saves none
#'   of them, keeping the production per-run files small.
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

  arms <- list()

  # -- nuisances (random forest) ----------------------------------------------
  cat("Nuisances: RF double crossfit...\n")
  nz_double <- timed(nuisance_double_rf(X, Y, W, fold_indices, fold_pairs, num.threads))
  cat("Nuisances: RF single crossfit...\n")
  nz_single <- timed(nuisance_single_rf(X, Y, W, fold_indices, num.threads))
  cat("Nuisances: RF out-of-bag...\n")
  nz_oob <- timed(nuisance_oob_rf(X, Y, W, num.threads))

  # -- DR learner, random forest ---------------------------------------------
  cat("DR-RF variants...\n")

  s <- timed(stage2_crossfit_rf(X, nz_double$value$po, X_test, fold_indices, num.threads))
  arms$dcf <- arm("dr_rf", "dcf", s$value$tau, s$value$tau_test, nz_double$time, s$time,
                  s$value$tau_test_folds, truth_test)

  s <- timed(stage2_crossfit_rf(X, nz_single$value$po, X_test, fold_indices, num.threads))
  arms$scf_scf <- arm("dr_rf", "scf_scf", s$value$tau, s$value$tau_test, nz_single$time, s$time,
                      s$value$tau_test_folds, truth_test)

  s <- timed(stage2_whole_rf(X, nz_single$value$po, X_test, num.threads))
  arms$scf_oob <- arm("dr_rf", "scf_oob", s$value$tau_oob, s$value$tau_test,
                      nz_single$time, s$time, truth_test = truth_test,
                      variance = s$value$var_oob)

  s <- timed(stage2_whole_rf(X, nz_oob$value$po, X_test, num.threads))
  arms$oob_oob <- arm("dr_rf", "oob_oob", s$value$tau_oob, s$value$tau_test,
                      nz_oob$time, s$time, truth_test = truth_test,
                      variance = s$value$var_oob)

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
                          nz_single$time, s$time, truth_test = truth_test,
                          variance = s$value$var_oob)

  # grf's own internally cross-fit OOB nuisances - no separate nuisance stage
  s <- timed(cf_whole(X, Y, W, X_test, num.threads = num.threads))
  arms$cf_default <- arm("causal_forest", "cf_default", s$value$tau_oob, s$value$tau_test,
                         0, s$time, truth_test = truth_test,
                         variance = s$value$var_oob)
  # the hats grf cross-fit for itself, kept so cf_default's half-sample bootstrap
  # can hold them fixed like every other arm's (see the return value below)
  nz_cf_default <- list(Y.hat.cf = s$value$Y.hat.cf, W.hat = s$value$W.hat)

  # -- DR learner, SuperLearner ----------------------------------------------
  # no OOB analogue exists for SuperLearner, so the OOB arms are dropped from
  # this family
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
  }

  list(arms = arms, fold_indices = fold_indices,
       nuisances = list(nz_double = nz_double$value, nz_single = nz_single$value,
                        nz_oob = nz_oob$value, nz_cf_default = nz_cf_default))
}
