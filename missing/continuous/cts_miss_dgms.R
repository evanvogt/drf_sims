###############
# title: continuous DGMs with missingness and imputations
###############

# Define scenario parameters for continuous outcomes
continuous_scenario_params <- data.frame(
  scenario = 1:5,
  description = c(
    "No HTE",
    "Simple HTE - binary variable",
    "Two HTE variables",
    "Single effects + different interaction (X3 + X4 + X4*X5)",
    "Cosine HTE"
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

#' Introduce missingness to specified variables (continuous version)
introduce_missingness_continuous <- function(data, miss_vars, miss_prop, mech = "MAR", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(data)
  allvars <- colnames(data)
  indicators <- matrix(ncol = length(allvars), nrow = 0)
  for (miss in miss_vars) {
    ind <- as.numeric(!(allvars %in% miss))
    indicators <- rbind(indicators, ind)
  }
  if (length(miss_vars) > 1) {
    ind <- as.numeric(!(allvars %in% miss_vars))
    indicators <- rbind(indicators, ind)
  }

  data_with_miss <- ampute(data,
                           prop = miss_prop,
                           patterns = indicators)

  return(data_with_miss$amp)
}

#' Handle missing data using various approaches (continuous version)
handle_missing_data_continuous <- function(data, method = "complete_cases", m_imputations = 10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
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
               mean_val <- mean(data[[var]], na.rm = TRUE)
               imputed_data[[var]][is.na(data[[var]])] <- mean_val
               message(paste("Imputed", sum(is.na(data[[var]])), "missing values in", var, "with mean =", round(mean_val, 3)))
             } else if (is.logical(data[[var]]) || all(data[[var]] %in% c(0, 1, NA))) {
               mode_val <- as.numeric(names(sort(table(data[[var]]), decreasing = TRUE))[1])
               imputed_data[[var]][is.na(data[[var]])] <- mode_val
               message(paste("Imputed", sum(is.na(data[[var]])), "missing values in", var, "with mode =", mode_val))
             }
           }
           return(imputed_data)
         },
         "multiple_imputation" = {
           if (!requireNamespace("mice", quietly = TRUE)) {
             stop("Package 'mice' is required for multiple imputation. Please install it with: install.packages('mice')")
           }
           message(paste("Performing multiple imputation with", m_imputations, "imputations..."))
           mice_result <- suppressMessages(
             mice::mice(data, m = m_imputations, method = 'pmm', printFlag = FALSE, seed = seed)
           )
           imputed_datasets <- vector("list", m_imputations)
           for (i in 1:m_imputations) {
             imputed_datasets[[i]] <- mice::complete(mice_result, i)
           }
           final_data <- data
           missing_vars <- names(data)[sapply(data, function(x) any(is.na(x)))]
           for (var in missing_vars) {
             missing_indices <- which(is.na(data[[var]]))
             if (is.numeric(data[[var]])) {
               imputed_values <- sapply(missing_indices, function(idx) {
                 values <- sapply(imputed_datasets, function(imp_data) imp_data[[var]][idx])
                 median(values)
               })
             } else {
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

#' Generate continuous outcome dataset for specified scenario with optional missingness
generate_continuous_scenario_data <- function(scenario, n, return_truth = TRUE, seed = NULL,
                                              add_missingness = TRUE, miss_type = "prognostic", 
                                              miss_prop = 0.1, miss_seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (!scenario %in% 1:10) stop("Scenario must be between 1 and 5")
  if (add_missingness) {
    if (!miss_type %in% c("prognostic", "predictive", "both")) stop("miss_type must be 'prognostic', 'predictive', or 'both'")
    if (miss_prop < 0 || miss_prop > 1) stop("miss_prop must be between 0 and 1")
  }
  params <- continuous_scenario_params[continuous_scenario_params$scenario == scenario, ]
  s_total <- params$s_err + params$s2
  diff <- power.t.test(n = n/2, delta = NULL, sd = s_total, power = 0.75)$delta
  bW <- round(-diff, digits = 2)
  W <- rbinom(n, 1, 0.5)
  X1 <- rbinom(n, 1, params$X1_prob)
  X2 <- rnorm(n, 0, params$s2)
  X3 <- if(params$needs_X3) rbinom(n, 1, params$X3_prob) else NULL
  X4 <- if(params$needs_X4) rnorm(n, 0, params$s4) else NULL  
  X5 <- if(params$needs_X5) rnorm(n, 0, params$s5) else NULL
  err <- rnorm(n, 0, params$s_err)
  treatment_effect <- switch(scenario,
                             rep(bW, n),
                             bW + params$b3 * X3,
                             bW + params$b3 * X3 + params$b4 * X4,
                             bW + params$b3 * X3 + params$b4 * X4 + params$b45 * X4 * X5,
                             bW + params$b4 * cos(X4)
  )
  Y <- params$b0 + params$b1 * X1 + params$b2 * X2 + W * treatment_effect + err
  X01 <- rnorm(n, 0, 1)
  X02 <- rnorm(n, 0, 1) 
  X03 <- rnorm(n, 0, 1)
  cats <- sample(c("A", "B", "C"), size = n, replace = TRUE, prob = c(0.45, 0.3, 0.25))
  X04 <- as.integer(cats == "A")
  X05 <- as.integer(cats == "B")
  dataset_vars <- list(Y = Y, W = W, X1 = X1, X2 = X2)
  if(params$needs_X3) dataset_vars$X3 <- X3
  if(params$needs_X4) dataset_vars$X4 <- X4  
  if(params$needs_X5) dataset_vars$X5 <- X5
  dataset_vars <- c(dataset_vars, list(X01 = X01, X02 = X02, X03 = X03, X04 = X04, X05 = X05))
  dataset <- as.data.frame(dataset_vars)
  if (add_missingness) {
    prognostic_vars <- c("X2")
    predictive_vars <- c()
    if(params$needs_X3) predictive_vars <- c(predictive_vars, "X3")
    if(params$needs_X4) predictive_vars <- c(predictive_vars, "X4")
    if(params$needs_X5) predictive_vars <- c(predictive_vars, "X5")
    if(length(predictive_vars) != 0) predictive_vars <- sample(predictive_vars, 1)
    miss_vars <- switch(miss_type,
                        "prognostic" = prognostic_vars,
                        "predictive" = predictive_vars,
                        "both" = c(prognostic_vars, predictive_vars)
    )
    dataset <- introduce_missingness_continuous(dataset, miss_vars, miss_prop, miss_seed)
  }
  result <- list(dataset = dataset, bW = bW)
  if(return_truth) {
    p0 <- params$b0 + params$b1 * X1 + params$b2 * X2
    p1 <- params$b0 + params$b1 * X1 + params$b2 * X2 + treatment_effect
    tau <- p1 - p0
    truth <- data.frame(p0 = p0, p1 = p1, tau = tau)
    result$truth <- truth
  }
  return(result)
}

#' Comprehensive function to generate continuous data and handle missingness
generate_and_process_continuous_data <- function(scenario, n, return_truth = TRUE, seed = NULL,
                                                 add_missingness = FALSE, miss_type = "prognostic", 
                                                 miss_prop = 0.1, miss_seed = NULL,
                                                 handle_missing = "none", m_imputations = 10, 
                                                 imputation_seed = NULL) {
  data_result <- generate_continuous_scenario_data(
    scenario = scenario, n = n, return_truth = return_truth, seed = seed,
    add_missingness = add_missingness, miss_type = miss_type, 
    miss_prop = miss_prop, miss_seed = miss_seed
  )
  # Track original complete cases indices if required
  orig_missing_indices <- NULL
  if (handle_missing == "complete_cases" && add_missingness) {
    orig_missing_indices <- complete.cases(data_result$dataset)
  }
  # Handle missing data if requested
  if (handle_missing != "none" && add_missingness) {
    processed_dataset <- handle_missing_data_continuous(
      data_result$dataset, 
      method = handle_missing, 
      m_imputations = m_imputations,
      seed = imputation_seed
    )
    data_result$dataset <- processed_dataset
    data_result$missing_method <- handle_missing
    # Subset truth using the original complete case indices
    if (handle_missing == "complete_cases" && return_truth && !is.null(orig_missing_indices)) {
      data_result$truth <- data_result$truth[orig_missing_indices, , drop = FALSE]
    }
  }
  return(data_result)
}

#' Generate oracle formula and parameters for a continuous outcome scenario
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
