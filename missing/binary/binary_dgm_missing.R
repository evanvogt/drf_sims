###############
# title: Data Generating Functions - binary outcome with missingness
# author: Ellie Van Vogt (enhanced)
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

#' Introduce missingness to specified variables
#' 
#' @param data Dataset to introduce missingness to
#' @param miss_vars Character vector of variable names to make missing
#' @param miss_prop Proportion of values to make missing (0-1)
#' @param seed Optional seed for reproducibility
#' @return Dataset with missingness introduced
introduce_missingness <- function(data, miss_vars, miss_prop, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(data)
  data_with_miss <- data
  
  for (var in miss_vars) {
    if (var %in% colnames(data)) {
      # Randomly select indices to make missing
      miss_indices <- sample(1:n, size = floor(n * miss_prop))
      data_with_miss[[var]][miss_indices] <- NA
    } else {
      warning(paste("Variable", var, "not found in dataset"))
    }
  }
  
  return(data_with_miss)
}

#' Generate dataset for specified scenario with optional missingness
#' 
#' @param scenario Integer 1-10 specifying which scenario to generate
#' @param n Sample size
#' @param return_truth Logical, whether to return true values (p0, p1, tau)
#' @param seed Optional seed for reproducibility
#' @param add_missingness Logical, whether to add missingness
#' @param miss_type Character: "prognostic" (X1, X2), "predictive" (X3, X4, X5), or "both"
#' @param miss_prop Proportion of values to make missing (0-1)
#' @param miss_seed Optional separate seed for missingness (for reproducibility)
#' @return List with dataset and optionally truth values
generate_binary_scenario_data <- function(scenario, n, return_truth = TRUE, seed = NULL,
                                          add_missingness = FALSE, miss_type = "prognostic", 
                                          miss_prop = 0.1, miss_seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Check valid scenario
  if (!scenario %in% 1:10) {
    stop("Scenario must be between 1 and 10")
  }
  
  # Check valid missingness parameters
  if (add_missingness) {
    if (!miss_type %in% c("prognostic", "predictive", "both")) {
      stop("miss_type must be 'prognostic', 'predictive', or 'both'")
    }
    if (miss_prop < 0 || miss_prop > 1) {
      stop("miss_prop must be between 0 and 1")
    }
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
  
  # Add missingness if requested
  if (add_missingness) {
    # Define variable groups
    prognostic_vars <- c("X2") # missingness only in one variable
    predictive_vars <- c()
    
    # Add predictive variables that exist in this scenario
    if(params$needs_X3) predictive_vars <- c(predictive_vars, "X3")
    if(params$needs_X4) predictive_vars <- c(predictive_vars, "X4")
    if(params$needs_X5) predictive_vars <- c(predictive_vars, "X5")
    
    # Choose one of the prognostic variables (if selected) to be missing
    if(length(predictive_vars) != 0) predictive_vars <- sample(predictive_vars, 1)
    
    # Determine which variables to make missing
    miss_vars <- switch(miss_type,
                        "prognostic" = prognostic_vars,
                        "predictive" = predictive_vars,
                        "both" = c(prognostic_vars, predictive_vars)
    )
    
    # Introduce missingness
    dataset <- introduce_missingness(dataset, miss_vars, miss_prop, miss_seed)
  }
  
  result <- list(dataset = dataset, bW = bW)
  
  # True CATEs (computed on complete data before missingness)
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

#' Handle missing data using various approaches
#' 
#' @param data Dataset with missing data
#' @param method Character: "complete_cases", "mean_imputation", or "multiple_imputation"
#' @param m_imputations Number of imputations for multiple imputation (default: 10)
#' @param seed Optional seed for reproducibility
#' @return Dataset with missing data handled according to specified method
handle_missing_data <- function(data, method = "complete_cases", m_imputations = 10, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Check if there's any missing data
  if (!any(is.na(data))) {
    message("No missing data found. Returning original dataset.")
    return(data)
  }
  
  method <- match.arg(method, c("complete_cases", "mean_imputation", "multiple_imputation"))
  
  switch(method,
         "complete_cases" = {
           complete_data <- data[complete.cases(data), ]
           n_removed <- nrow(data) - nrow(complete_data)
           message(paste("Complete case analysis: Removed", n_removed, "observations with missing data"))
           message(paste("Final sample size:", nrow(complete_data)))
           return(complete_data)
         },
         
         "mean_imputation" = {
           imputed_data <- data
           missing_vars <- names(data)[sapply(data, function(x) any(is.na(x)))]
           
           for (var in missing_vars) {
             if (is.numeric(data[[var]])) {
               # Mean imputation for continuous variables
               mean_val <- mean(data[[var]], na.rm = TRUE)
               imputed_data[[var]][is.na(data[[var]])] <- mean_val
               message(paste("Imputed", sum(is.na(data[[var]])), "missing values in", var, "with mean =", round(mean_val, 3)))
             } else if (is.logical(data[[var]]) || all(data[[var]] %in% c(0, 1, NA))) {
               # Mode imputation for binary variables
               mode_val <- as.numeric(names(sort(table(data[[var]]), decreasing = TRUE))[1])
               imputed_data[[var]][is.na(data[[var]])] <- mode_val
               message(paste("Imputed", sum(is.na(data[[var]])), "missing values in", var, "with mode =", mode_val))
             }
           }
           return(imputed_data)
         },
         
         "multiple_imputation" = {
           # Check if mice package is available
           if (!requireNamespace("mice", quietly = TRUE)) {
             stop("Package 'mice' is required for multiple imputation. Please install it with: install.packages('mice')")
           }
           
           # Perform multiple imputation
           message(paste("Performing multiple imputation with", m_imputations, "imputations..."))
           
           # Suppress mice output for cleaner console
           mice_result <- suppressMessages(
             mice::mice(data, m = m_imputations, method = 'pmm', printFlag = FALSE, seed = seed)
           )
           
           # Extract all imputed datasets
           imputed_datasets <- vector("list", m_imputations)
           for (i in 1:m_imputations) {
             imputed_datasets[[i]] <- mice::complete(mice_result, i)
           }
           
           # For each variable with missing data, take median across imputations
           final_data <- data
           missing_vars <- names(data)[sapply(data, function(x) any(is.na(x)))]
           
           for (var in missing_vars) {
             missing_indices <- which(is.na(data[[var]]))
             
             if (is.numeric(data[[var]])) {
               # For continuous variables: median across imputations
               imputed_values <- sapply(missing_indices, function(idx) {
                 values <- sapply(imputed_datasets, function(imp_data) imp_data[[var]][idx])
                 median(values)
               })
             } else {
               # For binary variables: mode across imputations
               imputed_values <- sapply(missing_indices, function(idx) {
                 values <- sapply(imputed_datasets, function(imp_data) imp_data[[var]][idx])
                 as.numeric(names(sort(table(values), decreasing = TRUE))[1])
               })
             }
             
             final_data[[var]][missing_indices] <- imputed_values
             message(paste("Imputed", length(missing_indices), "missing values in", var, "using median/mode of", m_imputations, "imputations"))
           }
           
           return(final_data)
         }
  )
}

#' Comprehensive function to generate data and handle missingness
#' 
#' @param scenario Integer 1-10 specifying which scenario to generate
#' @param n Sample size
#' @param return_truth Logical, whether to return true values (p0, p1, tau)
#' @param seed Optional seed for reproducibility
#' @param add_missingness Logical, whether to add missingness
#' @param miss_type Character: "prognostic", "predictive", or "both"
#' @param miss_prop Proportion of values to make missing (0-1)
#' @param miss_seed Optional separate seed for missingness
#' @param handle_missing Character: "none", "complete_cases", "mean_imputation", or "multiple_imputation"
#' @param m_imputations Number of imputations for multiple imputation
#' @param imputation_seed Optional seed for imputation
#' @return List with processed dataset and optionally truth values
generate_and_process_data <- function(scenario, n, return_truth = TRUE, seed = NULL,
                                      add_missingness = FALSE, miss_type = "prognostic", 
                                      miss_prop = 0.1, miss_seed = NULL,
                                      handle_missing = "none", m_imputations = 10, 
                                      imputation_seed = NULL) {
  
  # Generate data with potential missingness
  data_result <- generate_binary_scenario_data(
    scenario = scenario, n = n, return_truth = return_truth, seed = seed,
    add_missingness = add_missingness, miss_type = miss_type, 
    miss_prop = miss_prop, miss_seed = miss_seed
  )
  
  # Handle missing data if requested
  if (handle_missing != "none" && add_missingness) {
    processed_dataset <- handle_missing_data(
      data_result$dataset, 
      method = handle_missing, 
      m_imputations = m_imputations,
      seed = imputation_seed
    )
    
    # Update the result
    data_result$dataset <- processed_dataset
    data_result$missing_method <- handle_missing
    
    # If using complete cases, need to subset truth values too
    if (handle_missing == "complete_cases" && return_truth) {
      original_n <- nrow(data_result$truth)
      final_n <- nrow(processed_dataset)
      if (final_n < original_n) {
        # Find which rows were kept (this assumes complete cases removes from the end, 
        # but mice keeps original row order)
        complete_indices <- complete.cases(
          generate_binary_scenario_data(scenario, n, FALSE, seed, add_missingness, miss_type, miss_prop, miss_seed)$dataset
        )
        data_result$truth <- data_result$truth[complete_indices, ]
      }
    }
  }
  
  return(data_result)
}