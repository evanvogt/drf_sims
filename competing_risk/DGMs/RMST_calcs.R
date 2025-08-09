find_detectable_rmst_diff <- function(n_per_arm = 150,
                                      target_power = 0.75,
                                      tau = 24,
                                      admin_censor_time = 36,
                                      alpha = 0.05,
                                      nsim = 500) {
  
  cat("Finding detectable RMST difference for n =", n_per_arm, "per arm\n")
  cat("Target power:", target_power, "\n")
  
  # Search range for RMST difference
  rmst_diffs <- seq(1, 8, by = 0.5)
  powers <- numeric(length(rmst_diffs))
  
  for (i in seq_along(rmst_diffs)) {
    cat("Testing RMST difference =", rmst_diffs[i], "months...")
    
    power_result <- rmst_power_simulation(
      n_per_arm = n_per_arm,
      rmst_diff = rmst_diffs[i],
      tau = tau,
      admin_censor_time = admin_censor_time,
      nsim = nsim,
      alpha = alpha
    )
    
    powers[i] <- power_result$power
    cat(" Power =", round(powers[i], 3), "\n")
    
    # Stop if we've exceeded target power
    if (powers[i] >= target_power) break
  }
  
  # Interpolate to find exact detectable difference
  if (any(powers >= target_power)) {
    # Linear interpolation
    idx <- which(powers >= target_power)[1]
    if (idx > 1) {
      # Interpolate between idx-1 and idx
      x1 <- rmst_diffs[idx - 1]
      x2 <- rmst_diffs[idx]
      y1 <- powers[idx - 1]
      y2 <- powers[idx]
      
      detectable_diff <- x1 + (target_power - y1) * (x2 - x1) / (y2 - y1)
    } else {
      detectable_diff <- rmst_diffs[idx]
    }
  } else {
    detectable_diff <- NA
    cat("Warning: Target power not achieved with tested range\n")
  }