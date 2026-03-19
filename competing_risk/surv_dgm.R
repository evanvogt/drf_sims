##########
# title: generating survival data using cumulative incidence (Beyersmenn 2009)
##########
require(future.apply)
require(dplyr)
# generating data from weibull distributions

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
  shape1_base = 1,
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

# helper functions
cum_haz_all <- function(t, shape1, shape2, scale1, scale2) {
  (t/scale1)^shape1 + (t/scale2)^shape2
}

find_time <- function(u, shape1, shape2, scale1, scale2) {
  f <- function(t) cum_haz_all(t, shape1, shape2, scale1, scale2) + log(u)
  uniroot(f, c(0, 200))$root
}

hazard <- function(t, shape, scale) {
  shape * (1/scale)^shape * t^(shape - 1)
}

truth_individual <- function(shape1, shape2, scale1_1, scale1_0, scale2_1, scale2_0, admin) {
  # hazards under each treatment
  h1_1 <- function(t) hazard(t, shape1, scale1_1)
  h1_0 <- function(t) hazard(t, shape1, scale1_0)
  
  h2_1 <- function(t) hazard(t, shape2, scale2_1)
  h2_0 <- function(t) hazard(t, shape2, scale2_0)
  
  # all-cause cumulative hazards
  cum_haz_all_1 <- function(t) cum_haz_all(t, shape1, shape2, scale1_1, scale2_1)
  cum_haz_all_0 <- function(t) cum_haz_all(t, shape1, shape2, scale1_0, scale2_0)
  
  # all cause survival functions
  S_all_1 <- function(t) exp(-cum_haz_all_1(t))
  S_all_0 <- function(t) exp(-cum_haz_all_0(t))
  
  # terms to integrate for getting CIFs
  
  f1_1 <- function(t) S_all_1(t) * h1_1(t)
  f1_0 <- function(t) S_all_0(t) * h1_0(t)
  
  f2_1 <- function(t) S_all_1(t) * h2_1(t)
  f2_0 <- function(t) S_all_0(t) * h2_0(t)
  
  # cause-specific cumulative incidence functions
  CIF1_1 <- Vectorize(function(t) integrate(f1_1, 0, t)$value)
  CIF1_0 <- Vectorize(function(t) integrate(f1_0, 0, t)$value)
  
  CIF2_1 <- Vectorize(function(t) integrate(f2_1, 0, t)$value)
  CIF2_0 <- Vectorize(function(t) integrate(f2_0, 0, t)$value)
  
  # cause-specific hazard estimands
  # RMTL - expected time from event i to the event horizon
  RMTL1_1 <- integrate(function(u) CIF1_1(u), 0, admin)$value
  RMTL1_0 <- integrate(function(u) CIF1_0(u), 0, admin)$value
  RMTL2_1 <- integrate(function(u) CIF2_1(u), 0, admin)$value
  RMTL2_0 <- integrate(function(u) CIF2_0(u), 0, admin)$value
  
  # RMSTc - expected time free from any event (essentially a composite outcome)
  RMSTc_1 <- integrate(function(u) 1 - CIF1_1(u) - CIF2_1(u), 0, admin)$value
  RMSTc_0 <- integrate(function(u) 1 - CIF1_0(u) - CIF2_0(u), 0, admin)$value
  
  # RMSTi - expected time to event i, irrespective of competing events 
  RMST1_1 <- admin - RMTL1_1
  RMST1_0 <- admin - RMTL1_0
  RMST2_1 <- admin - RMTL2_1
  RMST2_0 <- admin - RMTL2_0
  
  # subdistribution hazard estimands
  # sub-dist hazards
  sh_1 <- function(t) h1_1(t) / (1 + CIF2_1(t)/S_all_1(t))
  sh_0 <- function(t) h1_0(t) / (1 + CIF2_0(t)/S_all_0(t))
  
  # survival functions
  S_sh_1 <- Vectorize(function(t) exp(-integrate(sh_1, 0, t)$value))
  S_sh_0 <- Vectorize(function(t) exp(-integrate(sh_0, 0, t)$value))
  
  # RMST - time to event - keeping competing events in the risk set
  sh_RMST_1 <- integrate(function(u) S_sh_1(u), 0, admin)$value
  sh_RMST_0 <- integrate(function(u) S_sh_0(u), 0, admin)$value
  
  list(RMTL1_1 = RMTL1_1, RMTL1_0 = RMTL1_0, RMTL2_1 = RMTL2_1, RMTL2_0 = RMTL2_0,
       RMSTc_1 = RMSTc_1, RMSTc_0 = RMSTc_0, RMST1_1 = RMST1_1, RMST1_0 = RMST1_0,
       RMST2_1 = RMST2_1, RMST2_0 = RMST2_0, sh_RMST_1 = sh_RMST_1, sh_RMST_0 = sh_RMST_0)
}

generate_surv_data <- function(scenario, n, return_truth = TRUE, censoring = FALSE) {
  # scenario check
  if (!scenario %in% 1:7) {
    stop("Scenario must be between 1 and 7")
  }
  
  # Get parameters for this scenario
  params <- survival_scenario_params[survival_scenario_params$scenario == scenario, ]
  
  # Generate treatment and covariates
  W <- rbinom(n, 1, 0.5)
  X1 <- rbinom(n, 1, params$X1_prob)
  X2 <- rnorm(n, 0, 1)
  X3 <- rbinom(n, 1, params$X3_prob)
  
  # Prognostic effects on shape parameters (log-linear)
  log_shape1 <- log(params$shape1_base) + params$b1_shape1 * X1 + params$b2_shape1 * X2
  log_shape2 <- log(params$shape2_base) + params$b1_shape2 * X1 + params$b2_shape2 * X2
  
  shape1 <- exp(log_shape1)
  shape2 <- exp(log_shape2)
  
  # Treatment effects on scale parameters
  log_scale1 <- log(params$scale1_base) + W * (params$ATE1 * params$bW_scale1 + params$HTE1 * params$b3_scale1 * X3)
  
  log_scale2 <- log(params$scale2_base) + W * (params$ATE2 * params$bW_scale2 + params$HTE2 * params$b3_scale2 * X3)
  
  scale1 <- exp(log_scale1)
  scale2 <- exp(log_scale2)
  
  # Solve for survival time (joint hazard)
  u <- runif(n)
  Y <- future_mapply(find_time, u, shape1, shape2, scale1, scale2)
  
  # Generate cause
  h1 <- hazard(Y, shape1, scale1)
  h2 <- hazard(Y, shape2, scale2)
  
  prob <- h1/(h1 + h2)
  
  cause <- 2 - rbinom(n, 1, prob)
  
  # Administrative censoring (end of study)
  admin <- rep(params$event_horizon, n)
  D <- ifelse(Y > admin, 0, D)
  Y <- ifelse(Y > admin, admin, Y)
  
  # Uniform uniformative censoring (if required)
  if (censoring) {
    censor_time <- runif(n, 0, admin)
    D <- ifelse(Y > censor_time, 0, D)
    Y <- pmin(Y, censor_time)
  }
  
  # Add extra non-informative covariates
  X01 <- rnorm(n, 0, 1)
  X02 <- rnorm(n, 0, 1) 
  X03 <- rnorm(n, 0, 1)
  cats <- sample(c("A", "B", "C"), size = n, replace = TRUE, prob = c(0.45, 0.3, 0.25))
  X04 <- as.integer(cats == "A")
  X05 <- as.integer(cats == "B")
  
  # Dataset for analysis
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
  
  # truth calculation for RMST
  if (return_truth) {
    # scale under each treatment
    log_scale1_1 <- log(params$scale1_base) + (params$ATE1 * params$bW_scale1 + params$HTE1 * params$b3_scale1 * X3)
    
    scale1_1 <- exp(log_scale1_1)
    scale1_0 <- rep(params$scale1_base, n)
    
    log_scale2_1 <- log(params$scale2_base) + (params$ATE2 * params$bW_scale2 + params$HTE2 * params$b3_scale2 * X3)
    
    scale2_1 <- exp(log_scale2_1)
    scale2_0 <- rep(params$scale2_base, n)
    
    truth <- future_mapply(truth_individual, shape1, shape2, scale1_1, scale1_0, scale2_1, scale2_0, admin, SIMPLIFY = FALSE)
    truth <- do.call(rbind, truth) %>% data.frame()
    
    truth  <- truth %>%
      mutate(
        tau_RMTL1 = RMTL1_1 - RMTL1_0,
        tau_RMTL2 = RMTL2_1 - RMTL2_0,
        tau_RMSTc = RMSTc_1 - RMSTc_0,
        tau_RMST1 = RMST1_1 - RMST1_0,
        tau_RMST2 = RMST2_1 - RMST2_0,
        tau_sh = sh_RMST_1 - sh_RMST_0
      )
    result$truth <- truth
  }
  return(result)
}
