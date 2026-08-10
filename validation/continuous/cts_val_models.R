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
# The TE-VIM (treatment-effect variable importance) computation below is
# genuinely study-specific - nothing else in the repo uses it yet - so it
# stays here rather than moving into R/.

source(here::here("R", "cate_models.R"))

run_all_cate_methods <- function(data, n_folds = 10) {

  X <- as.matrix(data[, -c(1:2)])
  Y <- data$Y
  W <- data$W
  n_obs <- nrow(X)

  fold_indices <- sort(seq(n_obs) %% n_folds) + 1
  fold_list <- unique(fold_indices)
  fold_pairs <- utils::combn(fold_list, 2, simplify = FALSE)

  cat("Computing nuisance functions...\n")
  nuisances <- nuisance_rf(X, Y, W, fold_indices, fold_pairs)

  results <- list()

  cat("Running Causal Forest...\n")
  results$causal_forest <- run_causal_forest(X, Y, W, nuisances, fold_indices, fold_list)
  results$causal_forest$te_vims <- get_te_vims_causal_forest(
    X, Y, W, nuisances, results$causal_forest$tau, fold_indices, fold_list
  )

  cat("Running DR Random Forest...\n")
  results$dr_random_forest <- run_dr_random_forest(X, Y, W, nuisances, fold_indices, fold_list)
  results$dr_random_forest$te_vims <- get_te_vims(
    X, nuisances$po, results$dr_random_forest$tau, fold_indices, fold_list,
    nuisances$po_matrix
  )

  results
}

###################
# TE-VIMs (treatment-effect variable importance measures)
###################
# For each covariate, refit the second-stage regression on every fold with
# that covariate dropped, then compare its out-of-fold prediction error
# against the full model's - a pathwise-derivative-style importance measure.
# `get_te_vims` covers estimators whose second stage is a plain regression
# forest on a pseudo-outcome (DR Random Forest); `get_te_vims_causal_forest`
# refits a causal_forest instead, since that estimator's second stage takes
# (X, Y, W, nuisances) rather than (X, pseudo-outcome).

#' TE-VIMs for a DR-learner style estimator (pseudo-outcome regression)
#'
#' @param X covariate matrix
#' @param po pseudo-outcome vector (row means across folds)
#' @param tau full-model CATE estimates
#' @param fold_indices crossfitting fold assignment
#' @param fold_list unique fold ids
#' @param po_matrix optional n x V double-crossfit pseudo-outcome matrix;
#'   if NULL, `po` is used directly in every fold (as for the oracle arms)
get_te_vims <- function(X, po, tau, fold_indices, fold_list, po_matrix = NULL) {
  n_obs <- nrow(X)
  covariates <- colnames(X)

  sub_taus_list <- future_map(seq_along(covariates), function(i) {
    new_X <- as.matrix(X[, -i])
    covariate_sub_taus <- rep(NA, n_obs)

    for (fold in seq_along(fold_list)) {
      in_train <- fold_indices != fold
      in_fold <- !in_train

      DR_sub <- if (is.null(po_matrix)) {
        regression_forest(new_X[in_train, ], po[in_train])
      } else {
        regression_forest(new_X[in_train, ], po_matrix[in_train, fold])
      }

      covariate_sub_taus[in_fold] <- predict(DR_sub, newdata = new_X[in_fold, ])$predictions
    }

    covariate_sub_taus
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

#' TE-VIMs for the causal forest, refitting a causal_forest per dropped covariate
#'
#' @param nuisances output of R/cate_models.R::nuisance_rf() - Y.hat.cf_matrix
#'   and W.hat_matrix are the same nuisances run_causal_forest() itself uses
get_te_vims_causal_forest <- function(X, Y, W, nuisances, tau, fold_indices, fold_list) {
  n_obs <- nrow(X)
  covariates <- colnames(X)
  po <- nuisances$po

  sub_taus_list <- future_map(seq_along(covariates), function(i) {
    new_X <- as.matrix(X[, -i])
    covariate_sub_taus <- rep(NA, n_obs)

    for (fold in seq_along(fold_list)) {
      in_train <- fold_indices != fold
      in_fold <- !in_train

      forest <- causal_forest(new_X[in_train, ], Y[in_train], W[in_train],
                              nuisances$Y.hat.cf_matrix[in_train, fold],
                              nuisances$W.hat_matrix[in_train, fold])

      covariate_sub_taus[in_fold] <- predict(forest, newdata = new_X[in_fold, ])$predictions
    }

    covariate_sub_taus
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
