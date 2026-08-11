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
# The TE-VIM (treatment-effect variable importance) computation below is
# genuinely study-specific - nothing else in the repo uses it yet - so it
# stays here rather than moving into R/.

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

  cat("Running DR Random Forest...\n")
  results$dr_random_forest <- run_dr_random_forest(X, Y, W, nuisances)
  results$dr_random_forest$te_vims <- get_te_vims(
    X, nuisances$po, results$dr_random_forest$tau
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
