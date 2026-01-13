####################
# Title: simulating survival data from CSH model
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
  
  # Event horizon
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
  
  # Does treatment impact event?
  ATE1 = c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE),
  ATE2 = c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE),
  
  # Treatment effects on log-shape
  bW_shape1 = c(1, 0, 1, 1, 0, 1, 1),
  bW_shape2 = c(0, -1, 0, -1, -1, -1, -1),
  
  # Is there heterogeneity on the event?
  HTE1 = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE),
  HTE2 = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE),
  
  # predictive effects on log-shape
  # Event 1 shape effects (ICU discharge - positive values increase shape = faster discharge for treated)
  b3_shape1 = c(0, 0, 0.4, 0.4, 0, 0, 0.3),  # Only scenarios 2 and 4 have HTE on event 1
  
  # Event 2 shape effects (mortality - negative values decrease shape = lower early mortality risk)
  b3_shape2 = c(0, 0, 0, 0, -0.2, -0.15, -0.1), # Only scenarios 3 and 4 have HTE on event 2
  
  stringsAsFactors = FALSE
)

# useful functions for calculating survival quantities
hazard <- function(t, shape, scale) {
  shape * (1/scale)^shape * t^(shape - 1)
}

cum_hazard <- function(t, shape, scale) {
  (t/scale)^shape
}

generate_csh_data <- function(scenario, n, return_truth = TRUE, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Check valid scenario
  if (!scenario %in% 1:7) {
    stop("Scenario must be between 1 and 7")
  }
  
  # Get parameters for this scenario
  params <- survival_scenario_params[survival_scenario_params$scenario == scenario, ]
  
  # Generate basic variables 
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
  log_shape1 <- log(params$shape1_base) + W * (params$ATE1 * params$bW_shape1 + params$HTE1 * params$b3_shape1 * X3)
  
  log_shape2 <- log(params$shape2_base) + W * (params$ATE2 * params$bW_shape2 + params$HTE2 * params$b3_shape2 * X3)
  
  shape1 <- exp(log_shape1)
  shape2 <- exp(log_shape2)
  
  
  # generate survival time
  Y <- rep(NA, n)
  D <- rep(NA, n)
  for (i in 1:n) {
    # get all cause cumulative hazard
    lambda1 <- scale1[i]
    k1 <- shape1[i]
    lambda2 <- scale2[i]
    k2 <- shape2[i]
    
    cum_haz <- function(t) {
      (t/lambda1)^k1 + (t/lambda2)^k2
    }
    # inverse calculation for survival time
    u <- runif(1)
    
    f <- function(t)  {
      cum_haz(t) + log(u)
    }
    time <- uniroot(f, c(1e-8, 100))$root
    
    # generate the cause of the event
    h1 <- hazard(time, k1, lambda1)
    h2 <- hazard(time, k2, lambda2)
    
    prob <- h1/(h1 + h2)
    
    cause <- 2 - rbinom(1, 1, prob)
    
    # add to Y and D
    Y[i] <- time
    D[i] <- cause
  }
  
  # Administrative censoring
  admin <- params$event_horizon
  D <- ifelse(Y > admin, 0, D)
  Y <- ifelse(Y > admin, admin, Y)
  
  # Add uniform censoring?
  
  # Add irrelevant covariates
  X01 <- rnorm(n, 0, 1)
  X02 <- rnorm(n, 0, 1) 
  X03 <- rnorm(n, 0, 1)
  cats <- sample(c("A", "B", "C"), size = n, replace = TRUE, prob = c(0.45, 0.3, 0.25))
  X04 <- as.integer(cats == "A")
  X05 <- as.integer(cats == "B")
  
  # Compile dataset
  dataset <- data.frame(
    Y = Y,
    D = D,
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
  
  # calculate true cause-specific RMSTs
  if(return_truth) {
    # shape (affected by treatment)
    log_shape1_1 <- log(params$shape1_base) + (params$ATE1 * params$bW_shape1 + params$HTE1 * params$b3_shape1 * X3)
    
    shape1_1 <- exp(log_shape1_1)
    shape1_0 <- params$shape1_base
    
    log_shape2_1 <- log(params$shape2_base) + (params$ATE2 * params$bW_shape2 + params$HTE2 * params$b3_shape2 * X3)
    
    shape2_1 <- exp(log_shape2_1)
    shape2_0 <- params$shape2_base
    
    # scale (base and prognostic only)
    log_scale1 <- log(params$scale1_base) + params$b1_scale1 * X1 + params$b2_scale1 * X2
    log_scale2 <- log(params$scale2_base) + params$b1_scale2 * X1 + params$b2_scale2 * X2
    
    scale1 <- exp(log_scale1)
    scale2 <- exp(log_scale2)
    
    # rmst vectors to fill
    rmst_csh1_1 <- numeric(n)
    rmst_csh1_0 <- numeric(n)
    
    rmst_csh2_1 <- numeric(n)
    rmst_csh2_0 <- numeric(n)
    
    # intergrations for each participant
    for (i in 1:n) {
      # scale (lambda - not affected by treatment)
      lambda1 <- scale1[i]
      lambda2 <- scale2[i]
      
      # shape (k - affected by treatment)
      k1_1 <- shape1_1[i]
      k1_0 <- shape1_0
      k2_1 <- shape2_1[i]
      k2_0 <- shape2_0
      
      h1 <- function(t, W) {
        hazard(t, W*k1_1 + (1-W)*k1_0, lambda1)
      }
      h2 <- function(t, W) {
        hazard(t, W*k2_1 + (1-W)*k2_0, lambda2)
      }
      cum_haz_all <- function(t, W) {
        cum_hazard(t, W*k1_1 + (1-W)*k1_0, lambda1) + cum_hazard(t, W*k2_1 + (1-W)*k2_0, lambda2)
      }
      S_all <- function(t, W) {
        exp(-cum_haz_all(t, W))
      }
      
      # CS-CIF
      f1_1 <- function(t) {
        S_all(t, 1)*h1(t, 1)
      }
      f1_0 <- function(t) {
        S_all(t, 0)*h1(t, 0)
      }
      f2_1 <- function(t) {
        S_all(t, 1)*h2(t, 1)
      }
      f2_0 <- function(t) {
        S_all(t, 0)*h2(t, 0)
      }
      
      rmst_csh1_1[i] <- integrate(f1_1, 0, admin)$value
      rmst_csh1_0[i] <- integrate(f1_0, 0, admin)$value
      rmst_csh2_1[i] <- integrate(f2_1, 0, admin)$value
      rmst_csh2_0[i] <- integrate(f2_0, 0, admin)$value
    }
    
    tau_csh_1 <- rmst_csh1_1 - rmst_csh1_0
    tau_csh_2 <- rmst_csh2_1 - rmst_csh2_0
    
    truth <- data.frame(
      rmst_control_cs1 = rmst_csh1_0,
      rmst_treat_cs1   = rmst_csh1_1,
      tau_cs1          = tau_csh_1,
      
      rmst_control_cs2 = rmst_csh2_0,
      rmst_treat_cs2   = rmst_csh2_1,
      tau_cs2          = tau_csh_2,
      
      shape1_control = shape1_0,
      shape1_treat   = shape1_1,
      shape2_control = shape2_0,
      shape2_treat   = shape2_1
    )
    
    result$truth <- truth
  }
  
  return(result)
}
