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
  event_horizon = 30,
  
  # Baseline Weibull parameters for event 1
  shape1_base = 1.2,
  scale1_base = 10,
  
  # Baseline Weibull parameters for event 2
  shape2_base = 1,
  scale2_base = 45,
  
  # Prognostic effects on log-shape parameters
  b1_shape1 = -0.15,
  b2_shape1 = 0.25,
  b1_shape2 = -0.20,
  b2_shape2 = 0.40,
  
  # Does treatment impact event?
  ATE1 = c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE),
  ATE2 = c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE),
  
  # Treatment effects on log-scale
  bW_scale1 = c(-1, 0, -1, -1, 0, -1, -1),
  bW_scale2 = c(0, 0.7, 0, 0.7, 0.7, 0.7, 0.7),
  
  # Is there heterogeneity on the event?
  HTE1 = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE),
  HTE2 = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE),
  
  # predictive effects on log-scale
  b3_scale1 = c(0, 0, -0.6, -0.6, 0, 0, -0.5),
  b3_scale2 = c(0, 0, 0, 0, 0.4, 0.3, 0.25),
  
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
  W <- rbinom(n, 1, 0.5)
  X1 <- rbinom(n, 1, params$X1_prob)
  X2 <- rnorm(n, 0, 1)
  X3 <- rbinom(n, 1, params$X3_prob)
  
  # Calculate prognostic effects on shape parameters (log-linear)
  log_shape1 <- log(params$shape1_base) + params$b1_shape1 * X1 + params$b2_shape1 * X2
  log_shape2 <- log(params$shape2_base) + params$b1_shape2 * X1 + params$b2_shape2 * X2
  
  shape1 <- exp(log_shape1)
  shape2 <- exp(log_shape2)
  
  # Calculate treatment effects on scale parameters
  log_scale1 <- log(params$scale1_base) + W * (params$ATE1 * params$bW_scale1 + params$HTE1 * params$b3_scale1 * X3)
  
  log_scale2 <- log(params$scale2_base) + W * (params$ATE2 * params$bW_scale2 + params$HTE2 * params$b3_scale2 * X3)
  
  scale1 <- exp(log_scale1)
  scale2 <- exp(log_scale2)
  
  
  # generate survival time
  Y <- rep(NA, n)
  D <- rep(NA, n)
  for (i in 1:n) {
    # get all cause cumulative hazard
    k1 <- shape1[i]
    lambda1 <- scale1[i]
    k2 <- shape2[i]
    lambda2 <- scale2[i]
    
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
    # scale (affected by treatment)
    log_scale1_1 <- log(params$scale1_base) + (params$ATE1 * params$bW_scale1 + params$HTE1 * params$b3_scale1 * X3)
    
    scale1_1 <- exp(log_scale1_1)
    scale1_0 <- params$scale1_base
    
    log_scale2_1 <- log(params$scale2_base) + (params$ATE2 * params$bW_scale2 + params$HTE2 * params$b3_scale2 * X3)
    
    scale2_1 <- exp(log_scale2_1)
    scale2_0 <- params$scale2_base
    
    # shape (base and prognostic only)
    log_shape1 <- log(params$shape1_base) + params$b1_shape1 * X1 + params$b2_shape1 * X2
    log_shape2 <- log(params$shape2_base) + params$b1_shape2 * X1 + params$b2_shape2 * X2
    
    shape1 <- exp(log_shape1)
    shape2 <- exp(log_shape2)
    
    # rmst vectors to fill
    rmst_csh1_1 <- numeric(n)
    rmst_csh1_0 <- numeric(n)
    
    rmst_csh2_1 <- numeric(n)
    rmst_csh2_0 <- numeric(n)
    
    # integrations for each participant
    for (i in 1:n) {
      # shape (k - not affected by treatment)
      k1 <- shape1[i]
      k2 <- shape2[i]
      
      # scale (lambda - affected by treatment)
      lambda1_1 <- scale1_1[i]
      lambda1_0 <- scale1_0
      lambda2_1 <- scale2_1[i]
      lambda2_0 <- scale2_0
      
      # hazards for events under each treatment
      h1_1 <- function(t) hazard(t, k1, lambda1_1)
      h1_0 <- function(t) hazard(t, k1, lambda1_0)
      
      h2_1 <- function(t) hazard(t, k2, lambda2_1)
      h2_0 <- function(t) hazard(t, k2, lambda2_0)
      
      # all-cause cumulative hazards
      cum_haz_all_1 <- function(t) cum_hazard(t, k1, lambda1_1) + cum_hazard(t, k2, lambda2_1)
      cum_haz_all_0 <- function(t) cum_hazard(t, k1, lambda1_0) + cum_hazard(t, k2, lambda2_0)
      
      # all cause survival functions
      S_all_1 <- function(t) exp(-cum_haz_all_1(t))
      S_all_0 <- function(t) exp(-cum_haz_all_0(t))
      
      f1_1 <- function(t) S_all_1(t) * h1_1(t)
      f1_0 <- function(t) S_all_0(t) * h1_0(t)
      
      f2_1 <- function(t) S_all_1(t) * h2_1(t)
      f2_0 <- function(t) S_all_0(t) * h2_0(t)
      
      # cumulative incidence functions up to time t
      CIF1_1 <- function(t) integrate(f1_1, 0, t)$value
      CIF1_0 <- function(t) integrate(f1_0, 0, t)$value
      
      CIF2_1 <- function(t) integrate(f2_1, 0, t)$value
      CIF2_0 <- function(t) integrate(f2_0, 0, t)$value
      
      # RMST integration (Fubini's theorem to remove double intergral with CIFs)
      
      rmst_csh1_1[i] <- admin - integrate(function(u) (admin - u)*f1_1(u), 0, admin)$value
      rmst_csh1_0[i] <- admin - integrate(function(u) (admin - u)*f1_0(u), 0, admin)$value
      rmst_csh2_1[i] <- admin - integrate(function(u) (admin - u)*f2_1(u), 0, admin)$value
      rmst_csh2_0[i] <- admin - integrate(function(u) (admin - u)*f2_0(u), 0, admin)$value
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
      
      scale1_control = scale1_0,
      scale1_treat   = scale1_1,
      scale2_control = scale2_0,
      scale2_treat   = scale2_1
    )
    
    result$truth <- truth
  }
  
  return(result)
}