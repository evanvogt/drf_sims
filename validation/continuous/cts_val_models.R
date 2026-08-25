##########
# title: CATE estimation for the interim-analysis validation study
##########
# The nuisance estimation and the two estimators (causal_forest,
# dr_random_forest) live in R/cate_models.R - this study uses only those two,
# because its question is whether subgroups/variance/variable-importance found
# on one chunk of a trial replicate on the next, not which estimator is best.
#
# Adopting R/cate_models.R's run_causal_forest also switches this study's
# causal-forest arm to the shared Y.hat.cf nuisance (a regression of Y on X
# alone) rather than the old validation-local nuisance (Y on (W,X) with the
# observed W plugged back in). That is a deliberate change, not an accident -
# see validation/README.md's Status section.
#
# Both estimators (and the TE-VIM refits below) are whole-sample OOB now, not
# fold-crossfit - see R/cate_models.R's crossfitting-strategy note. n_folds is
# kept as a parameter only so the analysis driver's call site doesn't need to
# change; it is otherwise unused.
#
# The two variable-importance computations below (TE-VIM and surrogate
# TreeSHAP) and the subgroup interaction test are genuinely study-specific -
# nothing else in the repo uses them yet - so they stay here rather than moving
# into R/.

library(xgboost)
library(SHAPforxgboost)

source(here::here("R", "cate_models.R"))

run_all_cate_methods <- function(data, n_folds = 10) {

  X <- as.matrix(data[, -c(1:2)])
  Y <- data$Y
  W <- data$W

  cat("Computing nuisance functions...\n")
  nuisances <- nuisance_rf(X, Y, W)

  results <- list()

  cat("Running Causal Forest...\n")
  results$causal_forest <- run_causal_forest(X, Y, W, nuisances)
  results$causal_forest$te_vims <- get_te_vims_causal_forest(
    X, Y, W, nuisances$po, results$causal_forest$tau
  )
  results$causal_forest$shap_vims <- get_shap_vims(X, results$causal_forest$tau)

  cat("Running DR Random Forest...\n")
  results$dr_random_forest <- run_dr_random_forest(X, Y, W, nuisances)
  results$dr_random_forest$te_vims <- get_te_vims(
    X, nuisances$po, results$dr_random_forest$tau
  )
  results$dr_random_forest$shap_vims <- get_shap_vims(
    X, results$dr_random_forest$tau
  )

  results
}

###################
# TE-VIMs (treatment-effect variable importance measures)
###################
# For each covariate, refit the whole-sample stage-2 model with that
# covariate dropped and take its OOB predictions, then compare their
# prediction error against the full model's OOB predictions - a
# pathwise-derivative-style importance measure. Whole-sample OOB throughout,
# matching R/cate_models.R's stage2_whole_rf / run_causal_forest - this used
# to refit per fold instead, before those moved off fold-crossfitting.
# `get_te_vims` covers estimators whose second stage is a plain regression
# forest on a pseudo-outcome (DR Random Forest); `get_te_vims_causal_forest`
# refits a causal_forest instead, since that estimator takes (X, Y, W)
# directly rather than (X, pseudo-outcome).

#' TE-VIMs for a DR-learner style estimator (pseudo-outcome regression)
#'
#' @param X covariate matrix
#' @param po pseudo-outcome vector
#' @param tau full-model OOB CATE estimates
get_te_vims <- function(X, po, tau) {
  n_obs <- nrow(X)
  covariates <- colnames(X)

  sub_taus_list <- future_map(seq_along(covariates), function(i) {
    new_X <- as.matrix(X[, -i])
    predict(regression_forest(new_X, po))$predictions
  }, .options = furrr_options(seed = TRUE))

  sub_taus <- do.call(cbind, sub_taus_list)
  colnames(sub_taus) <- covariates

  r_tau <- (po - tau)^2

  te_vims <- apply(sub_taus, 2, function(sub_tau) {
    r_subtau <- (po - sub_tau)^2
    tevim <- sum(r_subtau - r_tau) / n_obs
    infl <- r_subtau - r_tau - tevim
    std_err <- sqrt(sum(infl^2)) / n_obs
    list(tevim = tevim, std_err = std_err)
  }) %>% simplify2array()

  as.data.frame(te_vims)
}

#' TE-VIMs for the causal forest, refitting a whole-sample causal_forest per
#' dropped covariate
#'
#' Each dropped-covariate forest cross-fits its own nuisances internally, same
#' as the full model (R/cate_models.R::run_causal_forest) - this used to
#' refit fold-wise against externally-supplied double-crossfit nuisances.
#'
#' @param po the po field of nuisance_rf()'s output, used for scoring only -
#'   the forest itself no longer takes externally-supplied nuisances
#' @param tau full-model OOB CATE estimates
get_te_vims_causal_forest <- function(X, Y, W, po, tau) {
  n_obs <- nrow(X)
  covariates <- colnames(X)

  sub_taus_list <- future_map(seq_along(covariates), function(i) {
    new_X <- as.matrix(X[, -i])
    predict(causal_forest(new_X, Y, W))$predictions
  }, .options = furrr_options(seed = TRUE))

  sub_taus <- do.call(cbind, sub_taus_list)
  colnames(sub_taus) <- covariates

  r_tau <- (po - tau)^2

  te_vims <- apply(sub_taus, 2, function(sub_tau) {
    r_subtau <- (po - sub_tau)^2
    tevim <- sum(r_subtau - r_tau) / n_obs
    infl <- r_subtau - r_tau - tevim
    std_err <- sqrt(sum(infl^2)) / n_obs
    list(tevim = tevim, std_err = std_err)
  }) %>% simplify2array()

  as.data.frame(te_vims)
}


###################
# Surrogate TreeSHAP
###################
# The second variable-importance measure. Where the TE-VIMs above ask "how much
# worse does the second stage predict without this covariate", TreeSHAP asks
# "how much of this unit's estimated CATE is attributable to this covariate",
# and averages |that| over units. The two need not agree, which is the point of
# carrying both: cts_val_analysis.R compares each one's chunk-1 vs chunk-2
# ranking separately.
#
# This is "Strategy 3" / indirect SHAP from Svensson et al.'s SHAP_CATE - fit a
# tuned xgboost surrogate to the *estimated* CATEs, then take exact TreeSHAP of
# the surrogate. Because it reads only (X, tau) it applies unchanged to both
# estimators, unlike the TE-VIMs which need an estimator-specific refit.
#
# Adapted from SHAP_CATE-main/FUNCTIONS/functions.R::cvboost3(), with four
# deliberate departures, all for runtime - this runs four times per array job
# across 1100 jobs:
#   - 3 learning rates x 3 depths (9 combos), not 8 x 3 (24)
#   - 5-fold CV and at most 2000 rounds, not 10-fold and 10000
#   - the tuning grid is parallelised over the workers the analysis driver has
#     already planned, with nthread = 1 inside each fit
#   - the xgb.cv prediction callback is dropped. Only the final xgb.train model
#     is needed for SHAP and the round count comes off the evaluation log, so
#     saving the fold models was pure memory - and it also avoids the callback
#     having been renamed cb.cv.predict -> xgb.cb.cv.predict at xgboost 2.0,
#     which would otherwise tie this to a version nothing in this repo pins.

CVBOOST_ETAS <- c(0.01, 0.05, 0.1)
CVBOOST_DEPTHS <- 2:4

#' Cross-validated xgboost surrogate for a continuous target
#'
#' @param x feature matrix
#' @param y target vector (here, the estimated CATEs)
#' @return the fitted xgb.Booster at the best (eta, max_depth, nrounds)
cvboost_cate <- function(x, y, k_folds = 5, ETAs = CVBOOST_ETAS,
                         TR.DEPTH = CVBOOST_DEPTHS, ntrees_max = 2000,
                         early_stopping_rounds = 50) {
  tune_grid <- expand.grid(eta = ETAs, max_depth = TR.DEPTH)

  base_param <- function(eta, max_depth) {
    list(objective = "reg:squarederror", eval_metric = "rmse",
         subsample = 0.9, colsample_bytree = 0.9,
         eta = eta, max_depth = max_depth, nthread = 1)
  }

  fits <- future_map(seq_len(nrow(tune_grid)), function(i) {
    param <- base_param(tune_grid$eta[i], tune_grid$max_depth[i])
    cvfit <- xgboost::xgb.cv(
      data = xgboost::xgb.DMatrix(data = x, label = y),
      params = param,
      nfold = k_folds,
      nrounds = ntrees_max,
      early_stopping_rounds = early_stopping_rounds,
      maximize = FALSE,
      verbose = FALSE
    )
    # which.min rather than $best_iteration: same answer, no dependence on a
    # field whose presence has moved around between xgboost versions
    losses <- as.data.frame(cvfit$evaluation_log)[["test_rmse_mean"]]
    list(param = param, loss = min(losses), nrounds = which.min(losses))
  }, .options = furrr_options(seed = TRUE))

  best <- fits[[which.min(vapply(fits, function(f) f$loss, numeric(1)))]]

  xgboost::xgb.train(
    data = xgboost::xgb.DMatrix(data = x, label = y),
    params = best$param,
    nrounds = best$nrounds
  )
}

#' Surrogate-TreeSHAP variable importance for a fitted CATE
#'
#' Shaped to match get_te_vims()'s output so cts_val_analysis.R can read either
#' with the same `[1, ]` indexing: one column per covariate, in `colnames(X)`
#' order. Only the mean-|SHAP| summary is returned, not the n x p per-unit SHAP
#' matrix - nothing downstream uses the latter and it would add ~200KB per run
#' to cts_val_all.RDS.
#'
#' @param X covariate matrix
#' @param tau estimated CATEs to explain
get_shap_vims <- function(X, tau) {
  surrogate <- cvboost_cate(X, tau)

  shap <- SHAPforxgboost::shap.values(xgb_model = surrogate, X_train = X)

  # shap.values() returns mean_shap_score sorted descending by importance, not
  # in X's column order. Reindexing is not cosmetic: cts_val_analysis.R ranks
  # this against te_vims by position, so leaving it sorted would silently
  # scramble every rank.
  out <- as.data.frame(as.list(shap$mean_shap_score[colnames(X)]))
  colnames(out) <- colnames(X)
  rownames(out) <- "mean_abs_shap"
  out
}


###################
# Subgroup interaction test
###################

#' Interaction p-value for `Y ~ W * v` in the later chunk
#'
#' Indexes the coefficient by name, not position. The positional form is what
#' made bottom_pval report the *intercept* p-value (`pvals_bottom[1]`) rather
#' than the W-by-subgroup interaction, and it is fragile for a second reason:
#' when `v` is constant, lm drops the interaction and the coefficient matrix has
#' fewer than four rows, so `[4]` reads whatever happens to be there.
#'
#' @param v subgroup indicator or covariate to interact with treatment
#' @return the W:v p-value, or NA_real_ if `v` carries no contrast (e.g. rpart
#'   predicted no bottom10 leaf into chunk 2)
interaction_pval <- function(Y, W, v) {
  d <- data.frame(Y = Y, W = W, v = v)
  if (length(unique(na.omit(d$v))) < 2) return(NA_real_)
  co <- summary(lm(Y ~ W * v, data = d))$coefficients
  if (!"W:v" %in% rownames(co)) return(NA_real_)
  unname(co["W:v", 4])
}
