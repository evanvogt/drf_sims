###############
# title: Data Generating Functions - binary outcome
# author: Ellie Van Vogt
# Date started: 
###############

# Define scenario parameters
binary_scenario_params <- data.frame(
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
  # Base parameters (same for all scenarios)
  X1_prob = 0.4,
  X3_prob = 0.7,
  b0 = -0.4,
  b1 = 0.5,
  b2 = 0.5,
  # Variable scenario parameters
  b3 = c(NA, -0.4, NA, -0.4, NA, 0.2, 0.2, 0.2, 0.2, 0.2),
  b4 = c(NA, NA, 0.2, 0.3, NA, 0.5, 0.5, 0.5, 0.5, -0.1),
  b5 = c(NA, NA, NA, NA, NA, NA, -0.5, -0.5, NA, NA),
  b34 = c(NA, NA, NA, NA, -0.5, -0.5, NA, NA, NA, NA),
  b45 = c(NA, NA, NA, NA, NA, NA, -0.5, -0.5, NA, NA),
  # Variables needed for each scenario
  needs_X3 = c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
  needs_X4 = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  needs_X5 = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
  stringsAsFactors = FALSE
)

#' Generate dataset for specified scenario
#' 
#' @param scenario Integer 1-10 specifying which scenario to generate
#' @param n Sample size
#' @param return_truth Logical, whether to return true values (p0, p1, tau)
#' @param seed Optional seed for reproducibility
#' @return List with dataset and optionally truth values
generate_binary_scenario_data <- function(scenario, n, return_truth = TRUE, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Check valid scenario
  if (!scenario %in% 1:10) {
    stop("Scenario must be between 1 and 10")
  }
  
  # Get parameters for this scenario
  params <- binary_scenario_params[binary_scenario_params$scenario == scenario, ]
  
  # Calculate bW for adequate power
  p1_base <- plogis(params$b0)
  p2 <- power.prop.test(n/2, p2 = p1_base, power = 0.75)$p1
  bW <- round(qlogis(p2) - params$b0, digits = 2)
  
  # Treatment and prognostic variables
  W <- rbinom(n, 1, 0.5)
  X1 <- rbinom(n, 1, params$X1_prob)
  X2 <- rnorm(n, 0, 1)
  
  # Predictive variables
  X3 <- if(params$needs_X3) rbinom(n, 1, params$X3_prob) else NULL
  X4 <- if(params$needs_X4) rnorm(n, 0, 1) else NULL  
  X5 <- if(params$needs_X5) rnorm(n, 0, 1) else NULL
  
  # select treatment effect function
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
                             bW + params$b4 * exp(X4)
  )
  
  # Linear predictor and outcome
  lp <- params$b0 + params$b1 * X1 + params$b2 * X2 + W * treatment_effect
  prob <- plogis(lp)
  Y <- rbinom(n, 1, prob)
  
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
  
  # True CATEs
  if(return_truth) {
    p0 <- plogis(params$b0 + params$b1 * X1 + params$b2 * X2)
    p1 <- plogis(params$b0 + params$b1 * X1 + params$b2 * X2 + treatment_effect)
    tau <- p1 - p0
    
    truth <- data.frame(p0 = p0, p1 = p1, tau = tau)
    result$truth <- truth
  }
  
  return(result)
}

#' Generate oracle formula and parameters for a scenario
#' 
#' @param scenario Integer 1-10 specifying which scenario
#' @param bW Treatment effect coefficient (from generate_scenario_data)
#' @return List with formula string and parameter values
get_binary_oracle_info <- function(scenario, bW) {
  
  params <- binary_scenario_params[binary_scenario_params$scenario == scenario, ]
  
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
                        "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b4*exp(X$X4))"
  )
  
  return(list(fmla = formula_str, params = param_list))
}