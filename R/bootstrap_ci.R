##########
# title: half-sample bootstrap confidence intervals for CATE estimates
##########
# One copy of what was duplicated across confidence_intervals/continuous/,
# confidence_intervals/binary/ and missing/ci_example/.
#
# Method: draw B half-samples stratified by fold, refit the second stage on each,
# and form the "half-sample root" tau_full - tau_half. The intervals are
# simultaneous over units: each root is standardised by its own bootstrap SD, the
# maximum over units is taken within each draw, and the (1 - alpha/2) quantile of
# that maximum sets a single critical value. That is why S_star is a scalar.
#
# No rescaling constant appears anywhere, and that is not an oversight: the half
# sample is nested in the full sample, so Cov(tau_n, tau_n/2) ~ Var(tau_n) and
# Var(tau_n - tau_n/2) ~ (2^2g - 1) Var(tau_n) for an estimator converging at
# rate g. At the parametric g = 1/2 that factor is exactly 1. A forest converging
# more slowly than that has roots which understate the target variance, which is
# one mechanism behind the sub-nominal coverage these studies keep measuring.
#
# Two families of bootstrap live here:
#   cf_half_boot / rf_half_boot          - per-fold crossfit stage 2
#   cf_oob_half_boot / rf_oob_half_boot  - whole-sample stage 2, OOB predictions
# All four share an argument order and count so a single table-driven loop can
# drive them (see crossfitting/confidence_intervals/cf_ci_analysis.R).

require(grf)
require(future)
require(furrr)

#' Draw a fold-stratified half sample
#'
#' Half of each fold is selected, so every fold stays represented and the
#' crossfitting structure survives the resampling.
half_sample <- function(fold_list, fold_membership, fold_sizes, n_obs) {
  keep <- rep(FALSE, n_obs)
  for (fold in fold_list) {
    fold_obs <- fold_membership[[fold]]
    n_half <- floor(fold_sizes[fold] / 2)
    keep[sample(fold_obs, n_half, replace = FALSE)] <- TRUE
  }
  keep
}

#' Turn a matrix of half-sample roots into a simultaneous interval
#'
#' @param draws n x B matrix of (tau - tau_half)
#' @param tau the point estimates the interval is centred on
#' @param na.rm whether to skip NA roots. FALSE (the default) is the behaviour
#'   every fold-structured caller wants, since those fill every cell. The OOB
#'   bootstraps' "out-of-half only" band deliberately masks the in-half cells,
#'   which leaves each unit ~B/2 roots and makes the column max a supremum over
#'   ~n/2 rather than n units - a smaller sup, so a smaller S_star. That is a
#'   property of the construction, not of this function.
simultaneous_band <- function(draws, tau, alpha, na.rm = FALSE) {
  lambda_hat <- apply(draws, 1, var, na.rm = na.rm)
  draws_norm <- abs(draws) / sqrt(lambda_hat)
  col_max <- apply(draws_norm, 2, max, na.rm = na.rm)
  S_star <- quantile(col_max, 1 - (alpha / 2))

  margin <- sqrt(lambda_hat) * S_star

  list(hb_lb = tau - margin, hb_ub = tau + margin, draws = draws)
}

#' Half-sample bootstrap for the causal forest
#'
#' @param nuisances the crossfit nuisance object. Double-crossfit nuisances
#'   (nuisance_double_rf) carry fold-indexed Y.hat.cf_matrix / W.hat_matrix
#'   columns, reused unchanged so only the forest is refit on each half
#'   sample. Single-crossfit nuisances (nuisance_single_rf) carry only the
#'   un-suffixed vector fields Y.hat.cf / W.hat - detected the same way
#'   cf_foldwise (crossfitting/cf_models.R) detects its own Y.hat/W.hat shape.
#' @param tau point estimates from the full sample
#' @param CI_boot number of bootstrap draws
#' @param CI_sf sample.fraction handed to the half-sample forests
cf_half_boot <- function(X, Y, W, nuisances, tau, CI_boot = 200, CI_sf = 0.5,
                         alpha = 0.05, fold_indices, fold_list) {

  n_obs <- nrow(X)
  fold_membership <- lapply(fold_list, function(fold) which(fold_indices == fold))
  fold_sizes <- lengths(fold_membership)

  # prefer the matrix fields when present (every existing caller has them);
  # single-crossfit nuisances have only the vector fields
  Y.hat <- if (!is.null(nuisances$Y.hat.cf_matrix)) nuisances$Y.hat.cf_matrix else nuisances$Y.hat.cf
  W.hat <- if (!is.null(nuisances$W.hat_matrix)) nuisances$W.hat_matrix else nuisances$W.hat
  single <- is.vector(Y.hat)
  stopifnot(is.vector(Y.hat) == is.vector(W.hat))

  draws <- future_map(seq_len(CI_boot), function(b) {
    half_samples <- half_sample(fold_list, fold_membership, fold_sizes, n_obs)

    tau_half_results <- lapply(fold_list, function(fold) {
      in_train <- half_samples & (fold_indices != fold)
      in_fold <- fold_indices == fold

      y_hat <- if (single) Y.hat[in_train] else Y.hat[in_train, fold]
      w_hat <- if (single) W.hat[in_train] else W.hat[in_train, fold]

      half_cf <- causal_forest(X[in_train, ], Y[in_train], W[in_train],
                               y_hat, w_hat, sample.fraction = CI_sf)
      predict(half_cf, newdata = X[in_fold, ])$predictions
    })

    tau_half <- rep(NA, nrow(X))
    for (i in fold_list) tau_half[fold_membership[[i]]] <- tau_half_results[[i]]

    tau - tau_half
  }, .options = furrr_options(seed = TRUE))

  simultaneous_band(do.call(cbind, draws), tau, alpha)
}

#' Half-sample bootstrap for the DR-learner second stage
#'
#' @param po pseudo-outcomes: either the n x V double-crossfitting matrix or a
#'   plain n-vector (the oracle arms)
rf_half_boot <- function(X, Y, W, po, tau, CI_boot = 200, CI_sf = 0.5,
                         alpha = 0.05, fold_indices, fold_list) {

  n_obs <- nrow(X)
  fold_membership <- lapply(fold_list, function(fold) which(fold_indices == fold))
  fold_sizes <- lengths(fold_membership)
  single <- is.vector(po)

  draws <- future_map(seq_len(CI_boot), function(b) {
    half_samples <- half_sample(fold_list, fold_membership, fold_sizes, n_obs)

    tau_half_results <- lapply(fold_list, function(fold) {
      in_train <- half_samples & (fold_indices != fold)
      in_fold <- fold_indices == fold

      y_train <- if (single) po[in_train] else po[in_train, fold]
      half_rf <- regression_forest(X[in_train, ], y_train, sample.fraction = CI_sf)
      predict(half_rf, newdata = X[in_fold, ])$predictions
    })

    tau_half <- rep(NA, nrow(X))
    for (i in fold_list) tau_half[fold_membership[[i]]] <- tau_half_results[[i]]

    tau - tau_half
  }, .options = furrr_options(seed = TRUE))

  simultaneous_band(do.call(cbind, draws), tau, alpha)
}

# ---- whole-sample (OOB) arms ------------------------------------------------
#
# The arms whose stage 2 is one forest over the whole sample, scored by its own
# out-of-bag predictions, have no per-fold structure for the two bootstraps above
# to refit against. One forest is refit per draw on an unstratified half sample
# instead, and tau_half is assembled by oob_half_predict.

#' Draw an unstratified half sample
#'
#' The whole-sample/OOB arms have no fold structure to preserve, so the half
#' sample is a plain floor(n/2) draw without replacement. Kept separate from
#' half_sample() rather than folded into it behind a NULL fold_list, so the two
#' resampling schemes stay individually readable.
oob_half_sample <- function(n_obs) {
  keep <- rep(FALSE, n_obs)
  keep[sample(n_obs, floor(n_obs / 2), replace = FALSE)] <- TRUE
  keep
}

#' Assemble tau_half for a whole-sample forest
#'
#' In-half rows take the half forest's own OOB predictions - the same functional
#' the point estimate is, since predict() with no newdata restricts each row to
#' its out-of-bag trees. Out-of-half rows were never seen by the forest at all,
#' so a plain newdata prediction is already honest for them.
#'
#' Neither branch is contaminated by the row's own outcome. The one asymmetry is
#' tree count: an OOB prediction averages the (1 - sample.fraction) share of
#' trees while a newdata prediction averages all of them, which at grf's default
#' num.trees is second-order Monte Carlo noise beside the statistical variance.
oob_half_predict <- function(forest, X, keep) {
  tau_half <- numeric(length(keep))
  tau_half[keep]  <- predict(forest)$predictions
  tau_half[!keep] <- predict(forest, newdata = X[!keep, , drop = FALSE])$predictions
  tau_half
}

#' Both OOB bands from one set of half-sample roots
#'
#' There are two defensible ways to score an OOB half-sample refit, and they come
#' off identical resamples and identical forest fits - so this is a paired
#' contrast, not two independent bootstraps, and masking costs nothing.
#'
#'   hb_*     every unit scored in every draw (oob_half_predict as-is)
#'   hb_out_* only rows the half forest never saw; in-half cells masked to NA,
#'            which leaves each unit ~B/2 roots and makes S_star a sup over ~n/2
#'            rather than n units. sup of |N(0,1)| grows like sqrt(2 log m), so
#'            at n = 500 that predicts a band roughly 6% narrower. Whether it
#'            shows up once the roots are correlated is a question for the study,
#'            not something to assume here.
#'
#' @param roots n x B matrix of (tau - tau_half)
#' @param kept n x B logical, TRUE where the row was in that draw's half sample
oob_bands <- function(roots, kept, tau, alpha) {
  roots_out <- roots
  roots_out[kept] <- NA_real_

  band_all <- simultaneous_band(roots, tau, alpha)
  band_out <- simultaneous_band(roots_out, tau, alpha, na.rm = TRUE)

  list(hb_lb = band_all$hb_lb, hb_ub = band_all$hb_ub,
       hb_out_lb = band_out$hb_lb, hb_out_ub = band_out$hb_ub,
       draws = roots)
}

#' Half-sample bootstrap for a whole-sample (OOB) DR-learner second stage
#'
#' The OOB counterpart of rf_half_boot. Pseudo-outcomes are held fixed and
#' sliced, exactly as the crossfit bootstraps hold their nuisances fixed - only
#' the final-stage forest is refit.
#'
#' @param po pseudo-outcomes, an n-vector (a whole-sample stage 2 has no fold
#'   columns to choose between)
#' @param fold_indices,fold_list accepted and ignored. A whole-sample arm has no
#'   folds; these are here so all four bootstraps in this file share an argument
#'   order and count, the same reason rf_half_boot accepts Y and W it never uses.
rf_oob_half_boot <- function(X, Y, W, po, tau, CI_boot = 200, CI_sf = 0.5,
                             alpha = 0.05, fold_indices = NULL, fold_list = NULL) {

  n_obs <- nrow(X)
  stopifnot(is.vector(po), length(po) == n_obs)

  res <- future_map(seq_len(CI_boot), function(b) {
    keep <- oob_half_sample(n_obs)
    half_rf <- regression_forest(X[keep, ], po[keep], sample.fraction = CI_sf)
    list(root = tau - oob_half_predict(half_rf, X, keep), keep = keep)
  }, .options = furrr_options(seed = TRUE))

  oob_bands(do.call(cbind, lapply(res, `[[`, "root")),
            do.call(cbind, lapply(res, `[[`, "keep")),
            tau, alpha)
}

#' Half-sample bootstrap for a whole-sample (OOB) causal forest
#'
#' @param nuisances an object carrying vector Y.hat.cf / W.hat. Both callers
#'   arrive in that shape: cf_full_oob passes its single-crossfit nuisances, and
#'   cf_default passes the Y.hat/W.hat grf itself cross-fit on the full sample
#'   (cf_whole returns them for exactly this purpose). Holding those fixed keeps
#'   cf_default under the same "only the final-stage forest is refit" contract as
#'   every other arm - letting grf re-cross-fit on each half sample would put
#'   nuisance variability into this one arm's band and nobody else's.
#' @param fold_indices,fold_list accepted and ignored, see rf_oob_half_boot
cf_oob_half_boot <- function(X, Y, W, nuisances, tau, CI_boot = 200, CI_sf = 0.5,
                             alpha = 0.05, fold_indices = NULL, fold_list = NULL) {

  n_obs <- nrow(X)
  Y.hat <- nuisances$Y.hat.cf
  W.hat <- nuisances$W.hat
  stopifnot(is.vector(Y.hat), is.vector(W.hat))

  res <- future_map(seq_len(CI_boot), function(b) {
    keep <- oob_half_sample(n_obs)
    half_cf <- causal_forest(X[keep, ], Y[keep], W[keep],
                             Y.hat[keep], W.hat[keep], sample.fraction = CI_sf)
    list(root = tau - oob_half_predict(half_cf, X, keep), keep = keep)
  }, .options = furrr_options(seed = TRUE))

  oob_bands(do.call(cbind, lapply(res, `[[`, "root")),
            do.call(cbind, lapply(res, `[[`, "keep")),
            tau, alpha)
}

# ---- pooling bootstrap intervals across multiple imputations ----------------

#' Combine half-sample bootstrap intervals over multiply-imputed datasets
#'
#' Used only by missing/ci_example, which is the study asking how to build a CATE
#' interval when the data were multiply imputed. Three strategies are computed
#' side by side so they can be compared:
#'   *_pooled  empirical quantiles of the bootstrap replicates stacked across
#'             imputations
#'   *_mib     Rubin-style: within + between variance, critical value averaged
#'             over the per-imputation maxima
#'   *_hybrid  one variance and one critical value from the stacked draws
#'
#' `alpha` is an argument here. In the original it was a FREE VARIABLE read from
#' the global environment, which happened to work only because the analysis
#' script defined `alpha` at top level - the function would have silently used
#' the wrong level, or failed, in any other caller.
#'
#' @param res_list one per-imputation result from cate_methods()
#' @param model which model's estimates to pool
#' @param alpha two-sided level
combine_mi_ci <- function(res_list, model, alpha = 0.05) {
  res <- list()

  tau_mat <- do.call(cbind, lapply(res_list, function(x) x[[model]][["tau"]]))
  tau <- rowMeans(tau_mat)
  res$tau <- tau

  # causal-forest variance estimates, where the model produced them
  var_mat <- do.call(cbind, lapply(res_list, function(x) x[[model]][["variance"]]))
  if (!is.null(var_mat)) {
    res$cf_variance <- rowMeans(var_mat) +
      (1 + 1 / length(res_list)) * apply(tau_mat, 1, var)
  }

  # --- pooled: stack the bootstrap estimates themselves and take quantiles
  tau_bs_mat <- do.call(cbind, lapply(res_list, function(x) {
    x[[model]][["tau"]] - x[[model]][["draws"]]
  }))
  res$lb_pooled <- apply(tau_bs_mat, 1, quantile, probs = (alpha / 2))
  res$ub_pooled <- apply(tau_bs_mat, 1, quantile, probs = 1 - (alpha / 2))

  # --- mib: Rubin's rules on the bootstrap variances
  mi_b_list <- lapply(res_list, function(x) {
    draws <- x[[model]][["draws"]]
    lambda_hat <- apply(draws, 1, var)
    draws_norm <- abs(draws) / sqrt(lambda_hat)
    col_max <- apply(draws_norm, 2, max)
    list(lambda_hat = lambda_hat, col_max = col_max,
         S_star = quantile(col_max, 1 - (alpha / 2)))
  })

  lambda_mat <- do.call(cbind, lapply(mi_b_list, `[[`, "lambda_hat"))
  tot_var <- rowMeans(lambda_mat) +
    (1 + 1 / length(res_list)) * apply(tau_mat, 1, var)

  col_max_vec <- do.call(cbind, lapply(mi_b_list, `[[`, "col_max"))
  S_star <- quantile(col_max_vec, 1 - (alpha / 2))

  margin <- sqrt(tot_var) * S_star
  res$lb_mib <- tau - margin
  res$ub_mib <- tau + margin

  # --- hybrid: one variance and one critical value from the stacked draws
  draws_mat <- do.call(cbind, lapply(res_list, function(x) x[[model]][["draws"]]))
  lambda_hat <- apply(draws_mat, 1, var)
  draws_norm <- abs(draws_mat) / sqrt(lambda_hat)
  col_max <- apply(draws_norm, 2, max)
  S_star <- quantile(col_max, 1 - (alpha / 2))

  # NOTE: sqrt(lambda_hat * S_star), not sqrt(lambda_hat) * S_star as in the
  # other two strategies and in simultaneous_band(). Preserved as written, but
  # it looks like a typo - it takes the square root of the critical value.
  margin <- sqrt(lambda_hat * S_star)
  res$lb_hybrid <- tau - margin
  res$ub_hybrid <- tau + margin

  res
}

# ---- sample.fraction calibration --------------------------------------------

#' Choose the sample.fraction that gets the bootstrap closest to nominal coverage
#'
#' Outcome-agnostic despite having lived in a file called cts_ci_sf_calibration.R;
#' the binary calibration study already sourced it.
#'
#' Coverage is measured against `tau.hat` as a plug-in truth: pseudo-outcome
#' residuals are resampled within each fold column, tau is re-estimated, and the
#' interval is asked whether it covers the original estimate.
find_optimal_sf <- function(X, Y, W, nuisances_rf, tau.hat, fold_indices,
                            sf_grid = seq(0.05, 0.5, 0.05), n_sim = 50,
                            CI_boot = 100, alpha = 0.05, verbose = TRUE) {

  fold_list <- unique(fold_indices)
  n_obs <- length(fold_indices)
  target <- 1 - alpha

  # residuals around tau.hat; the NA structure of po_matrix is preserved
  po_residuals <- sweep(nuisances_rf$po_matrix, 1, tau.hat, FUN = "-")

  if (verbose) cat("Calibrating sample.fraction across", length(sf_grid),
                   "values,", n_sim, "simulations each...\n")

  # sequential outer loop - stage_2_rf and rf_half_boot already use the workers
  results_by_sf <- lapply(sf_grid, function(sf) {

    if (verbose) cat(" sf =", sf, "\n")

    sim_results <- lapply(seq_len(n_sim), function(b) {

      po_sim <- matrix(NA_real_, nrow = n_obs, ncol = length(fold_list))
      for (j in fold_list) {
        valid_rows <- which(fold_indices != j)
        valid_resid <- po_residuals[valid_rows, j]
        valid_resid <- valid_resid[is.finite(valid_resid)]
        po_sim[valid_rows, j] <- tau.hat[valid_rows] +
          sample(valid_resid, size = length(valid_rows), replace = TRUE)
      }

      tau.hat_sim <- stage_2_rf(X, po_sim, fold_indices, fold_list)

      boot_res <- rf_half_boot(
        X = X, Y = Y, W = W, po = po_sim, tau = tau.hat_sim,
        CI_boot = CI_boot, CI_sf = sf, alpha = alpha,
        fold_indices = fold_indices, fold_list = fold_list
      )

      list(coverage = mean(tau.hat >= boot_res$hb_lb & tau.hat <= boot_res$hb_ub),
           ci_width = mean(boot_res$hb_ub - boot_res$hb_lb))
    })

    coverage_vec <- sapply(sim_results, `[[`, "coverage")
    width_vec <- sapply(sim_results, `[[`, "ci_width")

    list(sf = sf,
         mean_coverage = mean(coverage_vec),
         sd_coverage = sd(coverage_vec),
         mean_ci_width = mean(width_vec))
  })

  coverage_curve <- dplyr::bind_rows(lapply(results_by_sf, as.data.frame))

  optimal_idx <- which.min(abs(coverage_curve$mean_coverage - target))
  optimal_sf <- coverage_curve$sf[optimal_idx]

  if (verbose) {
    cat("Optimal sf:", optimal_sf,
        "(mean coverage:", round(coverage_curve$mean_coverage[optimal_idx], 3), ")\n")
  }

  list(optimal_sf = optimal_sf, coverage_curve = coverage_curve, n_sim = n_sim)
}
