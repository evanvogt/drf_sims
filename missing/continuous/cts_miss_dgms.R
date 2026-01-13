###############
# title: continuous DGMs with missingness and imputations
###############
require(dplyr)
require(mice)
require(missForest)
require(VIM)

# Define scenario parameters
continuous_scenario_params <- data.frame(
  scenario = 1:5,
  description = c(
    "No HTE",
    "Simple HTE - binary variable (X3)",
    "Two HTE variables (X3 + X4)",
    "Single effects + different interaction (X3 + X4 + X4*X5)",
    "Non-linear HTE (cos(X4))"
  ),
  # Base parameters
  X1_prob = 0.4,
  X3_prob = 0.7,
  # Baseline and prognostic effects
  b0 = c(0.4, 0.2, 0.4, 1, 0.4),
  b1 = -0.05,
  b2 = c(2, 2, 2, 2, 1),
  # Variable scenario parameters for treatment effects
  b3 = c(NA, 2,  0.3, 2, NA),
  b4 = c(NA, NA, -1,  0.5, 1),
  b5 = c(NA, NA,  NA, -0.5, NA),
  b45 = c(NA, NA, NA, -0.5, NA),
  # Variance parameters
  s2 = 1,
  s4 = 1,
  s5 = 1,
  s_err = 0.5,
  # Variables needed for each scenario
  needs_X3 = c(FALSE, TRUE, TRUE, TRUE, FALSE),
  needs_X4 = c(FALSE, FALSE, TRUE, TRUE, TRUE),
  needs_X5 = c(FALSE, FALSE, FALSE, TRUE, FALSE),
  stringsAsFactors = FALSE
)

#' Generate continuous outcome dataset for specified scenario (missing data reduced scenarios)
#' 
#' @param scenario Integer 1-5 specifying which scenario to generate
#' @param n Sample size
#' @param return_truth Logical, whether to return true values (p0, p1, tau)
#' @param seed Optional seed for reproducibility
generate_continuous_scenario_data <- function(scenario, n, return_truth = TRUE, seed = NULL) {
  # Checks
  if (!is.null(seed)) set.seed(seed)
  if (!scenario %in% 1:5) { stop("Scenario must be between 1 and 5") }
  
  # Get parameters, total variance, and appropriate bW
  params <- continuous_scenario_params[continuous_scenario_params$scenario == scenario, ]
  s_total <- params$s_err + params$s2
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
                             bW + params$b3 * X3 + params$b4 * X4,
                             bW + params$b3 * X3 + params$b4 * X4 + params$b45 * X4 * X5,
                             bW + params$b4 * cos(X4)
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

#' Introduce missingness into the simulated dataset
#'
#' @param data Simulated dataset to be amputed
#' @param type Where the missingness should be introduced ("predictive", "prognostic", "both")
#' @param prop Proportion of missingness in (0,1)
#' @param mech Mechanism of missingness, either MAR (based on observed variables) or AUX (based on unseen auxillary variable)
#' @param seed Optional seed
introduce_missingness_continuous <- function(data, type, prop, mech, seed = NULL) {
  # Checks
  if (!is.null(seed)) set.seed(seed)
  if (!type %in% c("prognostic", "predictive", "both")) stop("type must be 'prognostic', 'predictive', or 'both'")
  if (prop < 0 || prop > 1) stop("miss_prop must be between 0 and 1")
  
  # Get sample size and covariates
  n <- nrow(data)
  orig <- colnames(data)
  keep <- data %>% select(all_of(c("Y", "W", "X01", "X02", "X03", "X04", "X05")))
  data <- data %>% select(-all_of(c("Y", "W", "X01", "X02", "X03", "X04", "X05")))
  covs <- colnames(data)
  prog_vars <- c("X1", "X2")
  pred_vars <- setdiff(covs, prog_vars)
  
  # Identify variables for amputation
  miss_vars <- switch(type,
                      "prognostic" = prog_vars,
                      "predictive" = pred_vars,
                      "both" = c(pred_vars, prog_vars))
  
  if (mech == "AUX") {
    U <- rnorm(n, 0, 1)
    covs <- c(covs, "U")
    data <- cbind(data, U)
  }
  
  # Generate missingness pattern
  if (length(miss_vars) > 1) {
    indicators <- expand.grid(rep(list(c(0,1)), length(miss_vars)))
    indicators <- indicators[!apply(indicators, 1, function(x) all(x == 1)), ] # remove row with all 1s
    colnames(indicators) <- miss_vars
    if(length(miss_vars) < length(covs)) {
      observed <- setdiff(covs, miss_vars)
      indicators[observed] <- 1
    }
    indicators <- indicators %>% select(all_of(covs))
    indicators <- indicators[!apply(indicators, 1, function(x) all(x == 0)), ] # remove row with all 0s
  }
  
  if (length(miss_vars) == 1) {
    indicators <- ifelse(covs == miss_vars, 1, 0)
    names(indicators) <- covs
  }
  

  
  # Variables that can influence the missingness (only for when there is an auxillary variable)
  if (mech == "AUX") {
    weights <- matrix(0, ncol = length(covs), nrow = if (is.null(nrow(indicators))) 1 else nrow(indicators))
    weights[,length(covs)] <- 1
  }
  
  # Amputation step
  result <- ampute(data,
                   prop = prop,
                   patterns = indicators,
                   weights = if (mech == "AUX") weights else NULL)
  
  # Put dataset back together
  data <- cbind(result$amp, keep)
  data <- data %>% select(all_of(orig))
  
  return(data)
}

#' Handle missing data in the simulated dataset
#' 
#' @param data Simulated dataset with missing data
#' @param method Missing data handling method (complete cases, mean imputation, missForsest imputation, regression imputation, missing indicator, IPW, or none)
#' @param seed Optional seed
handle_missingness_continuous <- function(data, method, seed = NULL) {
  # Checks
  if (!is.null(seed)) set.seed(seed)
  if (!any(is.na(data))) {
    message("No missing data found. Returning original dataset.")
    return(data)
  }
  method <- match.arg(method, c("complete_cases", "mean_imputation", "missforest", "regression", "missing_indicator", "IPW", "none"))
  
  switch(method,
         "complete_cases" = {
           complete_data <- data[complete.cases(data), ]
           retained_indices <- complete.cases(data)
           n_removed <- nrow(data) - nrow(complete_data)
           message(paste("Complete case analysis: Removed", n_removed))
           message(paste("Final sample size:", nrow(complete_data)))
           list(data = complete_data, retained_indices = retained_indices)
         },
         "mean_imputation" = {
           imputed_data <- apply(data, 2, function(x) {
             replace(x, is.na(x), mean(x, na.rm = TRUE))
           }) %>% as.data.frame()
           message("Mean imputation complete")
           list(data = imputed_data)
         },
         "missforest" = {
           keep <- data %>% select(all_of(c("Y", "W")))
           df <- data %>% select(-all_of(c("Y", "W")))
           df <- df %>%
             mutate(across(everything(), function(x) {
               unique_vals <- unique(x[!is.na(x)])
               bin <- length(unique_vals) == 2 && all(unique_vals %in% c(0, 1))
               if (bin) {factor(x, levels = c(0, 1))} else {x}
             }))
           mf_imputed <- missForest(df)
           imputed_df <- mf_imputed$ximp
           imputed_df <- as.data.frame(lapply(imputed_df, function(x) {
             if (is.factor(x)) {as.numeric(as.character(x))} else {x}
           }))
           imputed_data <- cbind(keep, imputed_df)
           message("Imputation with missForest complete")
           list(data = imputed_data)
         },
         "regression" = {
           keep <- data %>% select(all_of(c("Y", "W")))
           df <- data %>% select(-all_of(c("Y", "W")))
           
           miss <- names(df)[sapply(df, function(x) {any(is.na(x))})]
           complete <- setdiff(names(df), miss)
           imputed_df <- df
           for (var in miss) {
             fmla <- as.formula(paste(var, "~", paste(complete, collapse = " + ")))
             temp_imputed <- regressionImp(fmla, df)
             imputed_df[[var]] <- temp_imputed[[var]]
           }
           message("Imputation via regression complete")
           imputed_data <- cbind(keep, imputed_df)
           list(data = imputed_data)
         },
         "missing_indicator" = {
           miss <- names(data)[sapply(data, function(x) {any(is.na(x))})]
           imputed_data <- data
           for (var in miss) {
             ind_name <- paste0(var, "_missing")
             imputed_data[[ind_name]] <- ifelse(is.na(imputed_data[[var]]), 1, 0)
             imputed_data[[var]] <- ifelse(is.na(imputed_data[[var]]), 
                                           mean(imputed_data[[var]], na.rm = TRUE), 
                                           imputed_data[[var]])
           }
           message("Missing indicators + mean imputation completed")
           list(data = imputed_data)
         },
         "IPW" = {
           df <- data %>% select(-all_of(c("Y", "W")))
           
           miss <- names(df)[sapply(df, function(x) {any(is.na(x))})]
           complete <- setdiff(names(df), miss)
           df$cc <- ifelse(complete.cases(df), 1, 0)
           
           fmla <- as.formula(paste("cc ~", paste(complete, collapse = " + ")))
           miss_lr <- glm(fmla, family = binomial, data = df)
           
           retained_indices <- complete.cases(df)
           complete_data <- data[retained_indices, ]
           n_removed <- nrow(data) - nrow(complete_data)
           
           ipw <- 1 / miss_lr$fitted.values[retained_indices]
           
           message(paste("IPW: removed", n_removed, "observations"))
           message(paste("Final sample size:", nrow(complete_data)))
           list(data = complete_data, ipw = ipw, retained_indices = retained_indices)
         },
         "none" = {
           message("No missing data handling applied, original data set with missingness returned")
           list(data = data)
         })
}

#' Simulate data, add missingness, and handle missingness in one
#' 
#' @param scenario Integer 1-5 specifying which scenario to generate
#' @param n Sample size
#' @param return_truth Logical, whether to return true values (p0, p1, tau)
#' @param type Where the missingness should be introduced ("predictive", "prognostic", "both")
#' @param prop Proportion of missingness in (0,1)
#' @param mech Mechanism of missingness, either MAR (based on observed variables) or AUX (based on unseen auxillary variable)
#' @param method Missing data handling method (complete cases, mean imputation, missForsest imputation, regression imputation, missing indicator, IPW, or none)
#' @param seed Optional seed
generate_and_process_continuous_data <- function(scenario, n, return_truth = TRUE, type, prop, mech, method, seed = NULL) {
  # Data and truth generation step
  data_result <- generate_continuous_scenario_data(scenario, n, return_truth, seed)
  
  # Add missingness
  miss_dataset <- introduce_missingness_continuous(data_result$dataset, type, prop, mech, seed)
  
  # Handle missing data
  processed_dataset <- handle_missingness_continuous(miss_dataset, method, seed)
  data_result$dataset <- processed_dataset$data
  data_result$missing_method <- method
  
  # Indices removal if required
  if (method %in% c("complete_cases", "IPW") & return_truth) {
    retained_indices <- processed_dataset$retained_indices
    data_result$truth <- data_result$truth[retained_indices,]
  }
  
  # IPW inclusion if required
  if (method == "IPW") {
    data_result$ipw <- processed_dataset$ipw
  }
  return(data_result)
}

#' Generate oracle formula and parameters for a continuous outcome scenario
#' 
#' @param scenario Integer 1-5 specifying scenario
#' @param bW ATE value from generated data
get_continuous_oracle_info <- function(scenario, bW) {
  params <- continuous_scenario_params[continuous_scenario_params$scenario == scenario, ]
  param_list <- list(
    b0 = params$b0,
    b1 = params$b1, 
    b2 = params$b2,
    bW = bW
  )
  if(!is.na(params$b3)) param_list$b3 <- params$b3
  if(!is.na(params$b4)) param_list$b4 <- params$b4
  if(!is.na(params$b5)) param_list$b5 <- params$b5
  if(!is.na(params$b45)) param_list$b45 <- params$b45
  formula_str <- switch(scenario,
                        "b0 + b1*X$X1 + b2*X$X2 + W*bW",
                        "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3)", 
                        "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3 + b4*X$X4)",
                        "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3 + b4*X$X4 + b45*X$X4*X$X5)",
                        "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b4*cos(X$X4))"
  )
  return(list(fmla = formula_str, params = param_list))
}