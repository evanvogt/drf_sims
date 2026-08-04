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
simultaneous_band <- function(draws, tau, alpha) {
  lambda_hat <- apply(draws, 1, var)
  draws_norm <- abs(draws) / sqrt(lambda_hat)
  col_max <- apply(draws_norm, 2, max)
  S_star <- quantile(col_max, 1 - (alpha / 2))

  margin <- sqrt(lambda_hat) * S_star

  list(hb_lb = tau - margin, hb_ub = tau + margin, draws = draws)
}

#' Half-sample bootstrap for the causal forest
#'
#' @param nuisances the double-crossfit nuisance object; the fold-indexed
#'   Y.hat.cf_matrix / W.hat_matrix columns are reused unchanged, so only the
#'   forest is refit on each half sample
#' @param tau point estimates from the full sample
#' @param CI_boot number of bootstrap draws
#' @param CI_sf sample.fraction handed to the half-sample forests
cf_half_boot <- function(X, Y, W, nuisances, tau, CI_boot = 200, CI_sf = 0.5,
                         alpha = 0.05, fold_indices, fold_list) {

  n_obs <- nrow(X)
  fold_membership <- lapply(fold_list, function(fold) which(fold_indices == fold))
  fold_sizes <- lengths(fold_membership)

  draws <- future_map(seq_len(CI_boot), function(b) {
    half_samples <- half_sample(fold_list, fold_membership, fold_sizes, n_obs)

    tau_half_results <- lapply(fold_list, function(fold) {
      in_train <- half_samples & (fold_indices != fold)
      in_fold <- fold_indices == fold

      half_cf <- causal_forest(X[in_train, ], Y[in_train], W[in_train],
                               nuisances$Y.hat.cf_matrix[in_train, fold],
                               nuisances$W.hat_matrix[in_train, fold],
                               sample.fraction = CI_sf)
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
