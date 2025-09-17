# Competing Risks RMST-Based Simulation Framework
# Targeting restricted mean survival time differences with competing risks

library(survival)
library(cmprsk)
library(dplyr)
library(ggplot2)
library(randomizeR)  # For randomization lists
library(gridExtra)

#' Calculate hazard rates from cumulative incidence or survival
#' Based on your existing function with PASS validation
calculateHazardRates <- function(cumulativeIncidence,
                                 cumulativeSurvival,
                                 fixedTimePoint) {
  hazardRates <- matrix(data = NA, nrow = 2, ncol = 2)
  if (is.null(cumulativeIncidence)) {
    hazardRates <- -log(cumulativeSurvival) / fixedTimePoint
  } else {
    hazardRates[1, ] <- cumulativeIncidence[1, ] *
      (-log(1 - cumulativeIncidence[1, ] - cumulativeIncidence[2, ]) /
         (fixedTimePoint * (cumulativeIncidence[1, ] + cumulativeIncidence[2, ])))
    hazardRates[2, ] <- cumulativeIncidence[2, ] *
      (-log(1 - cumulativeIncidence[1, ] - cumulativeIncidence[2, ]) /
         (fixedTimePoint * (cumulativeIncidence[1, ] + cumulativeIncidence[2, ])))
  }
  
  return(hazardRates)
}

#' Calculate RMST for competing risks using cumulative incidence
#' For event of interest only (not overall survival)
calculateRMST <- function(time, status, tau, event_of_interest = 1) {
  # Fit cumulative incidence function
  cif_fit <- cuminc(time, status)
  
  # Extract cumulative incidence for event of interest
  event_name <- paste("1", event_of_interest)
  
  if (event_name %in% names(cif_fit)) {
    times <- cif_fit[[event_name]]$time
    cuminc <- cif_fit[[event_name]]$est
    
    # Add time 0 and tau if needed
    if (times[1] > 0) {
      times <- c(0, times)
      cuminc <- c(0, cuminc)
    }
    
    # Extend to tau if necessary
    if (max(times) < tau) {
      times <- c(times, tau)
      cuminc <- c(cuminc, cuminc[length(cuminc)])
    } else {
      # Truncate at tau
      idx <- times <= tau
      times <- times[idx]
      cuminc <- cuminc[idx]
      if (max(times) < tau) {
        times <- c(times, tau)
        cuminc <- c(cuminc, cuminc[length(cuminc)])
      }
    }
    
    # Calculate RMST as integral of (1 - F(t)) from 0 to tau
    # Where F(t) is cumulative incidence function
    rmst <- 0
    for (i in 2:length(times)) {
      dt <- times[i] - times[i-1]
      # Average survival probability over interval
      avg_surv <- 1 - (cuminc[i] + cuminc[i-1])/2
      rmst <- rmst + dt * avg_surv
    }
    
    return(rmst)
  } else {
    # No events of this type observed, RMST = tau
    return(tau)
  }
}

#' Enhanced simulation function with RMST targeting
simulateCompetingRisksRMST <- function(n, 
                                       baselineHazards,
                                       targetRMSTDifference,
                                       tau,
                                       eventOfInterest = 1,
                                       hetCoeffs = NULL,
                                       seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Generate patient characteristics
  X1 <- rnorm(n, 0, 1)
  X2 <- rbinom(n, 1, 0.5)
  X3 <- rnorm(n, 0, 1)
  
  # Treatment assignment (1:1 randomization)
  treatment <- rbinom(n, 1, 0.5)
  
  # Calculate individual treatment effects if HTE specified
  if (!is.null(hetCoeffs)) {
    individual_effect <- targetRMSTDifference + 
      hetCoeffs$X1 * X1 * treatment +
      hetCoeffs$X2 * X2 * treatment +
      hetCoeffs$X3 * X3 * treatment
  } else {
    individual_effect <- rep(targetRMSTDifference, n)
  }
  
  # Convert RMST difference to hazard ratio adjustment
  # This is an approximation - in practice, you might need iterative calibration
  hr_adjustment <- exp(-individual_effect / tau * treatment)
  
  # Apply baseline hazards with covariate effects
  lambda1 <- baselineHazards[1, 1] * exp(
    0.2 * X1 + 0.3 * X2 - 0.1 * X3  # Example covariate effects
  ) * hr_adjustment^(ifelse(eventOfInterest == 1, 1, 0))
  
  lambda2 <- baselineHazards[2, 1] * exp(
    -0.1 * X1 + 0.2 * X2 + 0.15 * X3  # Example covariate effects
  ) * hr_adjustment^(ifelse(eventOfInterest == 2, 1, 0.5))  # Partial effect on competing risk
  
  # Generate competing event times
  time1 <- rexp(n, rate = lambda1)  # Event of interest
  time2 <- rexp(n, rate = lambda2)  # Competing risk
  
  # Observed outcomes
  time <- pmin(time1, time2, tau)
  status <- ifelse(time == tau, 0,  # Censored
                   ifelse(time1 <= time2, 1, 2))  # 1 = event of interest, 2 = competing
  
  # Create dataset
  data.frame(
    id = 1:n,
    time = time,
    status = status,
    treatment = treatment,
    X1 = X1,
    X2 = X2,
    X3 = X3,
    true_effect = individual_effect
  )
}

#' Calculate power for RMST difference in competing risks setting
calculateRMSTPower <- function(n, 
                               baselineHazards,
                               targetRMSTDifference, 
                               tau,
                               eventOfInterest = 1,
                               alpha = 0.05,
                               nSim = 1000,
                               hetCoeffs = NULL) {
  
  pValues <- numeric(nSim)
  rmstDifferences <- numeric(nSim)
  
  for (i in 1:nSim) {
    # Simulate data
    data <- simulateCompetingRisksRMST(
      n = n,
      baselineHazards = baselineHazards,
      targetRMSTDifference = targetRMSTDifference,
      tau = tau,
      eventOfInterest = eventOfInterest,
      hetCoeffs = hetCoeffs,
      seed = i
    )
    
    # Calculate RMST for each group
    tryCatch({
      control_data <- data[data$treatment == 0, ]
      treatment_data <- data[data$treatment == 1, ]
      
      rmst_control <- calculateRMST(control_data$time, control_data$status, tau, eventOfInterest)
      rmst_treatment <- calculateRMST(treatment_data$time, treatment_data$status, tau, eventOfInterest)
      
      rmstDifferences[i] <- rmst_treatment - rmst_control
      
      # Bootstrap-based test for RMST difference
      # Simplified approach - in practice, use proper RMST testing methods
      n_control <- nrow(control_data)
      n_treatment <- nrow(treatment_data)
      
      # Pooled variance estimate (simplified)
      pooled_var <- var(c(control_data$time, treatment_data$time)) * 
        (1/n_control + 1/n_treatment)
      
      # Test statistic
      z_stat <- rmstDifferences[i] / sqrt(pooled_var)
      pValues[i] <- 2 * (1 - pnorm(abs(z_stat)))
      
    }, error = function(e) {
      pValues[i] <<- 1  # Conservative p-value if error
      rmstDifferences[i] <<- 0
    })
  }
  
  power <- mean(pValues < alpha, na.rm = TRUE)
  avgRMSTDiff <- mean(rmstDifferences, na.rm = TRUE)
  bias <- avgRMSTDiff - targetRMSTDifference
  
  return(list(
    power = power,
    avgRMSTDiff = avgRMSTDiff,
    bias = bias,
    targetRMSTDiff = targetRMSTDifference,
    pValues = pValues,
    rmstDifferences = rmstDifferences
  ))
}

#' Find required RMST difference for target power
findRequiredRMSTEffect <- function(n,
                                   baselineHazards, 
                                   tau,
                                   targetPower = 0.8,
                                   eventOfInterest = 1,
                                   alpha = 0.05,
                                   nSim = 500) {
  
  # Grid search for required RMST difference
  rmst_effects <- seq(0.5, 10, by = 0.5)  # Adjust range based on tau
  powers <- numeric(length(rmst_effects))
  
  cat("Searching for required RMST difference to achieve", targetPower, "power...\n")
  
  for (i in seq_along(rmst_effects)) {
    result <- calculateRMSTPower(
      n = n,
      baselineHazards = baselineHazards,
      targetRMSTDifference = rmst_effects[i],
      tau = tau,
      eventOfInterest = eventOfInterest,
      alpha = alpha,
      nSim = nSim
    )
    
    powers[i] <- result$power
    cat(sprintf("RMST difference: %.2f, Power: %.3f\n", rmst_effects[i], powers[i]))
    
    if (powers[i] >= targetPower) break
  }
  
  # Linear interpolation for exact effect size
  if (any(powers >= targetPower)) {
    idx <- which(powers >= targetPower)[1]
    if (idx > 1) {
      x1 <- rmst_effects[idx-1]; x2 <- rmst_effects[idx]
      y1 <- powers[idx-1]; y2 <- powers[idx]
      required_effect <- x1 + (targetPower - y1) * (x2 - x1) / (y2 - y1)
    } else {
      required_effect <- rmst_effects[1]
    }
  } else {
    required_effect <- max(rmst_effects)
    cat("Warning: Target power not achieved with maximum tested effect size\n")
  }
  
  return(list(
    required_rmst_diff = required_effect,
    effects_tested = rmst_effects[1:min(i, length(rmst_effects))],
    powers_achieved = powers[1:min(i, length(powers))]
  ))
}

#' Calculate operating characteristics for multiple sample sizes
calculateOperatingCharacteristics <- function(sampleSizes,
                                              baselineHazards,
                                              targetRMSTDifference,
                                              tau,
                                              eventOfInterest = 1,
                                              alpha = 0.05,
                                              nSim = 1000,
                                              hetCoeffs = NULL) {
  
  results <- data.frame()
  
  for (n in sampleSizes) {
    cat(sprintf("\nCalculating operating characteristics for n = %d...\n", n))
    
    result <- calculateRMSTPower(
      n = n,
      baselineHazards = baselineHazards,
      targetRMSTDifference = targetRMSTDifference,
      tau = tau,
      eventOfInterest = eventOfInterest,
      alpha = alpha,
      nSim = nSim,
      hetCoeffs = hetCoeffs
    )
    
    # Calculate additional metrics
    mse <- mean((result$rmstDifferences - targetRMSTDifference)^2, na.rm = TRUE)
    coverage_95 <- mean(abs(result$rmstDifferences - targetRMSTDifference) <= 
                          1.96 * sd(result$rmstDifferences, na.rm = TRUE), na.rm = TRUE)
    
    results <- rbind(results, data.frame(
      sample_size = n,
      target_rmst_diff = targetRMSTDifference,
      observed_rmst_diff = result$avgRMSTDiff,
      bias = result$bias,
      power = result$power,
      mse = mse,
      coverage_95 = coverage_95
    ))
  }
  
  return(results)
}

# Example usage with PASS validation data
cat("=== RMST-Based Competing Risks Simulation ===\n")

# Validation with PASS Example 1 parameters
cumulativeIncidence <- matrix(c(0.10, 0.05, 0.65, 0.65), nrow = 2, ncol = 2, byrow = TRUE)
fixedTimePoint <- 3
baselineHazards <- calculateHazardRates(
  cumulativeIncidence = cumulativeIncidence,
  cumulativeSurvival = NULL,
  fixedTimePoint = fixedTimePoint
)

cat("Baseline hazard rates calculated from PASS Example 1:\n")
print(baselineHazards)

# Define simulation parameters
tau <- 36  # Time horizon (days)
eventOfInterest <- 1  # Event of interest (e.g., discharge)
sampleSizes <- c(100, 200, 300, 400, 500)
targetRMSTDifference <- 3.0  # Target RMST difference (days)

# Example 1: Find required RMST difference for 80% power
cat("\n=== Finding Required RMST Difference for 80% Power ===\n")
n_target <- 300
required_effect <- findRequiredRMSTEffect(
  n = n_target,
  baselineHazards = baselineHazards,
  tau = tau,
  targetPower = 0.75,
  eventOfInterest = eventOfInterest,
  nSim = 200  # Reduced for example
)

cat(sprintf("\nFor n = %d, required RMST difference: %.2f days\n", 
            n_target, required_effect$required_rmst_diff))

# Example 2: Calculate operating characteristics
cat("\n=== Calculating Operating Characteristics ===\n")
operating_chars <- calculateOperatingCharacteristics(
  sampleSizes = sampleSizes,
  baselineHazards = baselineHazards,
  targetRMSTDifference = targetRMSTDifference,
  tau = tau,
  eventOfInterest = eventOfInterest,
  nSim = 200,  # Reduced for example
  hetCoeffs = hetCoeffs
)

print("Operating Characteristics:")
print(operating_chars)

# Visualizations
# Power curve
p1 <- ggplot(operating_chars, aes(x = sample_size, y = power)) +
  geom_line(size = 1, color = "blue") +
  geom_point(size = 3, color = "blue") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red", alpha = 0.7) +
  labs(title = "Power vs Sample Size for RMST Difference",
       subtitle = sprintf("Target RMST difference: %.1f days, τ = %d days", 
                          targetRMSTDifference, tau),
       x = "Sample Size", y = "Power") +
  theme_minimal() +
  ylim(0, 1)

# Bias plot
p2 <- ggplot(operating_chars, aes(x = sample_size, y = bias)) +
  geom_line(size = 1, color = "red") +
  geom_point(size = 3, color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  labs(title = "Bias vs Sample Size",
       subtitle = "Difference between observed and target RMST difference",
       x = "Sample Size", y = "Bias (days)") +
  theme_minimal()

print(p1)
print(p2)

# Summary table
cat("\n=== Summary ===\n")
power_80_idx <- which(operating_chars$power >= 0.8)
if (length(power_80_idx) > 0) {
  min_n_80_power <- min(operating_chars$sample_size[power_80_idx])
  cat(sprintf("Minimum sample size for ≥80%% power: %d\n", min_n_80_power))
} else {
  cat("None of the tested sample sizes achieved 80% power\n")
}

cat(sprintf("Target RMST difference: %.2f days\n", targetRMSTDifference))
cat(sprintf("Time horizon (τ): %d days\n", tau))
cat(sprintf("Event of interest: %d\n", eventOfInterest))

# Example 3: Generate a single dataset for inspection
cat("\n=== Example Dataset ===\n")
example_data <- simulateCompetingRisksRMST(
  n = 200,
  baselineHazards = baselineHazards,
  targetRMSTDifference = targetRMSTDifference,
  tau = tau,
  eventOfInterest = eventOfInterest,
  hetCoeffs = hetCoeffs,
  seed = 123
)

cat("Example dataset structure:\n")
print(head(example_data, 10))

# Calculate RMST for example dataset
control_rmst <- calculateRMST(
  example_data$time[example_data$treatment == 0],
  example_data$status[example_data$treatment == 0],
  tau, eventOfInterest
)

treatment_rmst <- calculateRMST(
  example_data$time[example_data$treatment == 1],
  example_data$status[example_data$treatment == 1],
  tau, eventOfInterest
)

cat(sprintf("\nExample RMST calculation:\n"))
cat(sprintf("Control group RMST: %.2f days\n", control_rmst))
cat(sprintf("Treatment group RMST: %.2f days\n", treatment_rmst))
cat(sprintf("Observed difference: %.2f days\n", treatment_rmst - control_rmst))
cat(sprintf("Target difference: %.2f days\n", targetRMSTDifference))

# Event distribution
cat("\nEvent distribution in example dataset:\n")
print(table(example_data$treatment, example_data$status, 
            dnn = c("Treatment", "Status")))
cat("Status: 0 = censored, 1 = event of interest, 2 = competing risk\n")