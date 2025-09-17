####################
# Title: simulating survival data
###################

# create dataframe of survival parameters
survival_scenario_params <- data.frame(
  scenario = 1:7,
  description = c(
    "ATE on EOI only",
    "ATE on CE only",
    "simple HTE on EOI, no ATE on CE",
    "Simple HTE on EOI, ATE on CE",
    "Simple HTE on CE, no ATE on EOI",
    "simple HTE on CE, ATE on EOI",
    "Simple HTE on both events"
  ),
  # Base parameters
  X1_prob = 0.4,
  X3_prob = 0.7,
  
  # Event horizon (time point for RMST calculation)
  event_horizon = 30, # days
  
  # Baseline Weibull parameters for event 1 (ICU discharge - positive outcome)
  shape1_base = 1.8,    # Shape > 1 gives increasing hazard (more likely over time)
  scale1_base = 12,     # days - median discharge around day 8-10
  
  # Baseline Weibull parameters for event 2 (mortality - negative outcome)  
  shape2_base = 0.8,    # Shape < 1 gives decreasing hazard (higher risk early)
  scale2_base = 45,     # days - mortality events spread over longer period
  
  # Prognostic effects on log-scale parameters
  b1_scale1 = -0.15,  # female effect on discharge (faster discharge, negative = shorter time)
  b2_scale1 = 0.25,   # severity score effect on discharge (higher severity = longer stay)
  b1_scale2 = -0.20,  # female effect on mortality (protective, negative = longer survival)
  b2_scale2 = 0.40,   # severity score effect on mortality (higher severity = higher risk)
  
  # predictive effects on log-shape
  # Event 1 shape effects (ICU discharge - positive values increase shape = faster discharge for treated)
  b3_shape1 = c(0, 0.4, 0, 0.3),  # Only scenarios 2 and 4 have HTE on event 1
  
  # Event 2 shape effects (mortality - negative values decrease shape = lower early mortality risk)
  b3_shape2 = c(0, 0, -0.2, -0.15), # Only scenarios 3 and 4 have HTE on event 2
  
  stringsAsFactors = FALSE
)

#' Calculate Restricted Mean Survival Time for Weibull distribution
#' 
#' @param shape Weibull shape parameter
#' @param scale Weibull scale parameter  
#' @param horizon Time horizon for RMST calculation
#' @return RMST value
weibull_rmst <- function(shape, scale, horizon) {
  # RMST = integral from 0 to horizon of S(t) dt
  # For Weibull: S(t) = exp(-(t/scale)^shape)
  
  integrate_survival <- function(t, shape, scale) {
    exp(-(t/scale)^shape)
  }
  
  integral_result <- integrate(integrate_survival, lower = 0, upper = horizon, 
                               shape = shape, scale = scale)
  
  return(integral_result$value)
}

generate_survival_scenario_data <- function(scenario, n, return_truth = TRUE, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Check valid scenario
  if (!scenario %in% 1:4) {
    stop("Scenario must be between 1 and 4")
  }
  
  # Get parameters for this scenario
  params <- survival_scenario_params[survival_scenario_params$scenario == scenario, ]
  
  # Generate basic variables (always present)
  W <- rbinom(n, 1, 0.5)           # Treatment assignment
  X1 <- rbinom(n, 1, params$X1_prob)  # Prognostic binary variable
  X2 <- rnorm(n, 0, 1)             # Prognostic continuous variable
  X3 <- rbinom(n, 1, params$X3_prob)  # Predictive variable
  
  # Calculate prognostic effects on scale parameters (log-linear)
  log_scale1 <- log(params$scale1_base) + params$b1_scale1 * X1 + params$b2_scale1 * X2
  log_scale2 <- log(params$scale2_base) + params$b1_scale2 * X1 + params$b2_scale2 * X2
  
  scale1 <- exp(log_scale1)
  scale2 <- exp(log_scale2)
  
  # Calculate treatment effects on shape parameters
  shape1_effect <- W * params$b3_shape1 * X3
  shape2_effect <- W * params$b3_shape2 * X3
  
  # Final shapes
  shape1 <- pmax(params$shape1_base + shape1_effect, 0.1) 
  shape2 <- pmax(params$shape2_base + shape2_effect, 0.1)
  
  # Generate competing event times
  U1 <- runif(n)
  U2 <- runif(n)
  time1 <- scale1 * (-log(U1))^(1/shape1) 
  time2 <- scale2 * (-log(U2))^(1/shape2) 
  
  # Determine first event
  time_obs <- pmin(time1, time2)
  event_type <- ifelse(time1 <= time2, 1, 2)
  
  # Apply administrative censoring
  admin_censor_time <- params$event_horizon
  censored <- time_obs > admin_censor_time
  time_obs[censored] <- admin_censor_time
  event_type[censored] <- 0
  
  # Add irrelevant covariates
  X01 <- rnorm(n, 0, 1)
  X02 <- rnorm(n, 0, 1) 
  X03 <- rnorm(n, 0, 1)
  cats <- sample(c("A", "B", "C"), size = n, replace = TRUE, prob = c(0.45, 0.3, 0.25))
  X04 <- as.integer(cats == "A")
  X05 <- as.integer(cats == "B")
  
  # Compile dataset
  dataset <- data.frame(
    Y = time_obs,
    D = event_type,
    W = W,
    X1 = X1,
    X2 = X2,
    X3 = X3,
    X01 = X01,
    X02 = X02,
    X03 = X03,
    X04 = X04,
    X05 = X05
  )
  
  result <- list(dataset = dataset)
  
  # --- Calculate truths if requested ---
  if(return_truth) {
    
    # Shapes under control/treatment (independent of subject’s W)
    shape1_control <- rep(params$shape1_base, n)
    shape2_control <- rep(params$shape2_base, n)
    
    shape1_treat <- pmax(params$shape1_base + params$b3_shape1 * X3, 0.1)
    shape2_treat <- pmax(params$shape2_base + params$b3_shape2 * X3, 0.1)
    
    rmst_control_subdist <- numeric(n)
    rmst_treat_subdist <- numeric(n)
    
    rmst_control_cs1 <- numeric(n) # Cause-specific (event 1 only)
    rmst_treat_cs1 <- numeric(n)
    
    rmst_control_cs2 <- numeric(n) # Cause-specific (event 2 only)
    rmst_treat_cs2 <- numeric(n)
    
    for(i in 1:n) {
      # --- Subdistribution (with competing risks) ---
      S1_c <- function(t) exp(-(t/scale1[i])^shape1_control[i])
      S2_c <- function(t) exp(-(t/scale2[i])^shape2_control[i])
      S_overall_c <- function(t) S1_c(t) * S2_c(t)
      
      rmst_control_subdist[i] <- integrate(S_overall_c, 0, params$event_horizon)$value
      
      S1_t <- function(t) exp(-(t/scale1[i])^shape1_treat[i])
      S2_t <- function(t) exp(-(t/scale2[i])^shape2_treat[i])
      S_overall_t <- function(t) S1_t(t) * S2_t(t)
      
      rmst_treat_subdist[i] <- integrate(S_overall_t, 0, params$event_horizon)$value
      
      # --- Cause-specific: event 1 only ---
      S1_c_only <- function(t) exp(-(t/scale1[i])^shape1_control[i])
      S1_t_only <- function(t) exp(-(t/scale1[i])^shape1_treat[i])
      rmst_control_cs1[i] <- integrate(S1_c_only, 0, params$event_horizon)$value
      rmst_treat_cs1[i]   <- integrate(S1_t_only, 0, params$event_horizon)$value
      
      # --- Cause-specific: event 2 only ---
      S2_c_only <- function(t) exp(-(t/scale2[i])^shape2_control[i])
      S2_t_only <- function(t) exp(-(t/scale2[i])^shape2_treat[i])
      rmst_control_cs2[i] <- integrate(S2_c_only, 0, params$event_horizon)$value
      rmst_treat_cs2[i]   <- integrate(S2_t_only, 0, params$event_horizon)$value
    }
    
    tau_subdist <- rmst_treat_subdist - rmst_control_subdist
    tau_cs1 <- rmst_treat_cs1 - rmst_control_cs1
    tau_cs2 <- rmst_treat_cs2 - rmst_control_cs2
    
    truth <- data.frame(
      # Subdistribution truth
      rmst_control_subdist = rmst_control_subdist,
      rmst_treat_subdist   = rmst_treat_subdist,
      tau_subdist          = tau_subdist,
      
      # Cause-specific truths
      rmst_control_cs1 = rmst_control_cs1,
      rmst_treat_cs1   = rmst_treat_cs1,
      tau_cs1          = tau_cs1,
      
      rmst_control_cs2 = rmst_control_cs2,
      rmst_treat_cs2   = rmst_treat_cs2,
      tau_cs2          = tau_cs2,
      
      # Shapes (just to track underlying params)
      shape1_control = shape1_control,
      shape1_treat   = shape1_treat,
      shape2_control = shape2_control,
      shape2_treat   = shape2_treat,
      
      # Predictive covariate
      X3 = X3
    )
    
    result$truth <- truth
  }
  
  return(result)
}
