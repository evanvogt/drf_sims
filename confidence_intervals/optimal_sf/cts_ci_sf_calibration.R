##########
# title: optimal sample.fraction calibration for half-sample bootstrap CIs
##########

require(grf)
require(dplyr)
require(future)
require(furrr)

find_optimal_sf <- function(
    X,
    Y,
    W,
    nuisances_rf,
    tau.hat,
    fold_indices,
    sf_grid   = seq(0.05, 0.5, 0.05),
    n_sim     = 50,
    CI_boot   = 100,
    alpha     = 0.05,
    verbose   = TRUE
) {

  fold_list <- unique(fold_indices)
  n_obs     <- length(fold_indices)
  target    <- 1 - alpha

  # Empirical residuals around tau.hat — NA structure of po_matrix is preserved
  po_residuals <- sweep(nuisances_rf$po_matrix, 1, tau.hat, FUN = "-")

  if (verbose) cat("Calibrating sample.fraction across", length(sf_grid),
                   "values,", n_sim, "simulations each...\n")

  # Sequential outer loop over sf values — inner stage_2_rf/rf_half_boot use future workers
  results_by_sf <- lapply(sf_grid, function(sf) {

    if (verbose) cat(" sf =", sf, "\n")

    sim_results <- lapply(seq_len(n_sim), function(b) {

      # Resample residuals within each fold column, preserving NA pattern
      po_sim <- matrix(NA_real_, nrow = n_obs, ncol = length(fold_list))
      for (j in fold_list) {
        valid_rows  <- which(fold_indices != j)
        valid_resid <- po_residuals[valid_rows, j]
        valid_resid <- valid_resid[is.finite(valid_resid)]
        po_sim[valid_rows, j] <- tau.hat[valid_rows] +
          sample(valid_resid, size = length(valid_rows), replace = TRUE)
      }

      # Re-estimate tau from simulated pseudo-outcomes
      tau.hat_sim <- stage_2_rf(X, po_sim, fold_indices, fold_list)

      # Bootstrap CIs with candidate sf
      boot_res <- rf_half_boot(
        X            = X,
        Y            = Y,
        W            = W,
        po           = po_sim,
        tau          = tau.hat_sim,
        CI_boot      = CI_boot,
        CI_sf        = sf,
        alpha        = alpha,
        fold_indices = fold_indices,
        fold_list    = fold_list
      )

      # Coverage against original tau.hat as plug-in truth
      coverage <- mean(tau.hat >= boot_res$hb_lb & tau.hat <= boot_res$hb_ub)
      ci_width <- mean(boot_res$hb_ub - boot_res$hb_lb)

      list(coverage = coverage, ci_width = ci_width)
    })

    coverage_vec <- sapply(sim_results, `[[`, "coverage")
    width_vec    <- sapply(sim_results, `[[`, "ci_width")

    list(
      sf            = sf,
      mean_coverage = mean(coverage_vec),
      sd_coverage   = sd(coverage_vec),
      mean_ci_width = mean(width_vec)
    )
  })

  coverage_curve <- bind_rows(lapply(results_by_sf, as.data.frame))

  optimal_idx <- which.min(abs(coverage_curve$mean_coverage - target))
  optimal_sf  <- coverage_curve$sf[optimal_idx]

  if (verbose) {
    cat("Optimal sf:", optimal_sf,
        "(mean coverage:", round(coverage_curve$mean_coverage[optimal_idx], 3), ")\n")
  }

  return(list(
    optimal_sf     = optimal_sf,
    coverage_curve = coverage_curve,
    n_sim          = n_sim
  ))
}
