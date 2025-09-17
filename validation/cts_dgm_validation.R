###############
# title: Streamlined data generating process for continuous outcomes multiple trials
###############

# Define scenario parameters for continuous outcomes
continuous_scenario_params <- data.frame(
  scenario = 1:10,
  description = c(
    "No HTE",
    "Simple HTE - binary variable",
    "Simple HTE - continuous variable", 
    "Two HTE variables",
    "Continuous-binary interaction (X3*X4)",
    "Single effects + interaction (X3 + X4 + X3*X4)",
    "Continuous-continuous interaction (X4*X5)",
    "Single effects + different interaction (X3 + X4 + X4*X5)",
    "Cosine HTE",
    "Exponential HTE"
  ),
  # Base parameters
  X1_prob = 0.4,
  X3_prob = 0.7,
  # Baseline and prognostic effects
  b0 = c(0.4, 0.2, 0.3, 0.4, 0.4, 1, 1, 1, 0.4, 0.4),
  b1 = -0.05,
  b2 = c(2, 2, 2, 2, 2, 2, 2, 2, 1, 2),
  # Variable scenario parameters for treatment effects
  b3 = c(NA, 2, NA, 0.3, NA, 2, 2, 2, NA, 0.3),
  b4 = c(NA, NA, -1, -1, NA, 0.5, 0.5, 0.5, 1, 0.1),
  b5 = c(NA, NA, NA, NA, NA, NA, -0.5, -0.5, NA, NA),
  b34 = c(NA, NA, NA, NA, 1, -0.5, NA, NA, NA, NA),
  b45 = c(NA, NA, NA, NA, NA, NA, -0.5, -0.5, NA, NA),
  # Variance parameters
  s2 = 1,
  s4 = 1,
  s5 = 1,
  s_err = 0.5,
  # Variables needed for each scenario
  needs_X3 = c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
  needs_X4 = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  needs_X5 = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
  stringsAsFactors = FALSE
)

#' Generate continuous outcome dataset for specified scenario
#' 
#' @param scenario Integer 1-10 specifying which scenario to generate
#' @param n Sample size
#' @param return_truth Logical, whether to return true values (p0, p1, tau)
#' @param seed Optional seed for reproducibility
#' @return List with dataset and optionally truth values
generate_continuous_scenario_data <- function(scenario, n, return_truth = TRUE, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Check valid scenario
  if (!scenario %in% 1:10) {
    stop("Scenario must be between 1 and 10")
  }
  
  # Get parameters for this scenario
  params <- continuous_scenario_params[continuous_scenario_params$scenario == scenario, ]
  
  # Total variance for power calculation
  s_total <- params$s_err + params$s2
  
  # Calculate bW for adequate power
  diff <- power.t.test(n = n/2, delta = NULL, sd = s_total, power = 0.75)$delta
  bW <- round(-diff, digits = 2)
  
  # Treatment and prognostic variables
  W <- rbinom(n, 1, 0.5)
  X1 <- rbinom(n, 1, params$X1_prob)
  X2 <- rnorm(n, 0, params$s2)
  
  # Predictive variables
  X3 <- if(params$needs_X3) rbinom(n, 1, params$X3_prob) else NULL
  X4 <- if(params$needs_X4) rnorm(n, 0, params$s4) else NULL  
  X5 <- if(params$needs_X5) rnorm(n, 0, params$s5) else NULL
  
  # Error term
  err <- rnorm(n, 0, params$s_err)
  
  # Select treatment effect function
  treatment_effect <- switch(scenario,
                             rep(bW, n),
                             bW + params$b3 * X3,
                             bW + params$b4 * X4,
                             bW + params$b3 * X3 + params$b4 * X4,
                             bW + params$b34 * X3 * X4,
                             bW + params$b3 * X3 + params$b4 * X4 + params$b34 * X3 * X4,
                             bW + params$b45 * X4 * X5,
                             bW + params$b3 * X3 + params$b4 * X4 + params$b45 * X4 * X5,
                             bW + params$b4 * cos(X4),
                             bW + params$b3 * X3 + params$b4 * exp(-abs(X4))
  )
  
  # Outcome calc
  Y <- params$b0 + params$b1 * X1 + params$b2 * X2 + W * treatment_effect + err
  
  # Unrelated variables
  X01 <- rnorm(n, 0, 1)
  X02 <- rnorm(n, 0, 1) 
  X03 <- rnorm(n, 0, 1)
  cats <- sample(c("A", "B", "C"), size = n, replace = TRUE, prob = c(0.45, 0.3, 0.25))
  X04 <- as.integer(cats == "A")
  X05 <- as.integer(cats == "B")
  
  # Build dataset 
  dataset_vars <- list(Y = Y, W = W, X1 = X1, X2 = X2)
  
  if(params$needs_X3) dataset_vars$X3 <- X3
  if(params$needs_X4) dataset_vars$X4 <- X4  
  if(params$needs_X5) dataset_vars$X5 <- X5
  
  # Extra variables
  dataset_vars <- c(dataset_vars, list(X01 = X01, X02 = X02, X03 = X03, X04 = X04, X05 = X05))
  
  dataset <- as.data.frame(dataset_vars)
  
  result <- list(dataset = dataset, bW = bW)
  
  # true CATEs
  if(return_truth) {
    p0 <- params$b0 + params$b1 * X1 + params$b2 * X2
    p1 <- params$b0 + params$b1 * X1 + params$b2 * X2 + treatment_effect
    tau <- p1 - p0
    
    truth <- data.frame(p0 = p0, p1 = p1, tau = tau)
    result$truth <- truth
  }
  
  return(result)
}

#' Generate oracle formula and parameters for a continuous outcome scenario
#' 
#' @param scenario Integer 1-10 specifying which scenario
#' @param bW Treatment effect coefficient (from generate_continuous_scenario_data)
#' @return List with formula string and parameter values
get_continuous_oracle_info <- function(scenario, bW) {
  
  params <- continuous_scenario_params[continuous_scenario_params$scenario == scenario, ]
  
  # Base parameters
  param_list <- list(
    b0 = params$b0,
    b1 = params$b1, 
    b2 = params$b2,
    bW = bW
  )
  
  # Scenario-specific parameters
  if(!is.na(params$b3)) param_list$b3 <- params$b3
  if(!is.na(params$b4)) param_list$b4 <- params$b4
  if(!is.na(params$b5)) param_list$b5 <- params$b5
  if(!is.na(params$b34)) param_list$b34 <- params$b34
  if(!is.na(params$b45)) param_list$b45 <- params$b45
  
  # Formula string
  formula_str <- switch(scenario,
                        "b0 + b1*X$X1 + b2*X$X2 + W*bW",
                        "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3)", 
                        "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b4*X$X4)",
                        "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3 + b4*X$X4)",
                        "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b34*X$X3*X$X4)",
                        "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3 + b4*X$X4 + b34*X$X3*X$X4)",
                        "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b45*X$X4*X$X5)",
                        "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3 + b4*X$X4 + b45*X$X4*X$X5)",
                        "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b4*cos(X$X4))",
                        "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3 + b4*exp(-abs(X$X4)))"
  )
  
  return(list(fmla = formula_str, params = param_list))
}