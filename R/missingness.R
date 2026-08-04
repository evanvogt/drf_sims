##########
# title: shared missingness generation and handling
##########
# One copy of what was duplicated across missing/continuous/, missing/binary/ and
# missing/ci_example/. The three differed only in the function-name suffix
# (_continuous / _binary), the MNAR-vs-AUX spelling, and the number of multiple
# imputations - none of which is about the outcome type, so the suffixes are gone
# and the imputation count is an argument.

require(dplyr)
require(mice)
require(missForest)
require(VIM)

source(here::here("R", "dgm_scenarios.R"))

# columns that are never amputated: the outcome, the treatment, and the
# deliberately-unrelated covariates
NEVER_MISSING <- c("Y", "W", "X01", "X02", "X03", "X04", "X05")
PROGNOSTIC_VARS <- c("X1", "X2")

MISS_METHODS <- c("complete_cases", "mean_imputation", "missforest", "regression",
                  "missing_indicator", "IPW", "multiple_imputation", "none")

#' Introduce missingness into a simulated dataset
#'
#' Builds an mice::ampute pattern matrix over the covariates being amputated and
#' applies it. Under MAR the missingness depends on the observed covariates;
#' under MNAR / MNAR-Y it is driven entirely by the unobserved U, which is what
#' the weight matrix below encodes.
#'
#' @param data simulated dataset
#' @param type which covariates to amputate: "prognostic", "predictive" or "both"
#' @param prop proportion of missingness, in (0, 1)
#' @param mech "MAR", "MNAR" or "MNAR-Y" ("AUX"/"AUX-Y" accepted as synonyms)
#' @param U the unobserved variable, required for the MNAR mechanisms
introduce_missingness <- function(data, type, prop, mech, U = NULL) {

  mech <- sub("^AUX", "MNAR", mech)

  if (!type %in% c("prognostic", "predictive", "both")) {
    stop("type must be 'prognostic', 'predictive', or 'both'")
  }
  if (prop < 0 || prop > 1) stop("miss_prop must be between 0 and 1")
  mnar <- mech %in% c("MNAR", "MNAR-Y")
  if (mnar && is.null(U)) {
    stop("unobserved variable U required for MNAR missingness generation")
  }

  orig <- colnames(data)
  keep <- data %>% select(all_of(NEVER_MISSING))
  data <- data %>% select(-all_of(NEVER_MISSING))
  covs <- colnames(data)

  prog_vars <- PROGNOSTIC_VARS
  pred_vars <- setdiff(covs, prog_vars)

  miss_vars <- switch(type,
                      "prognostic" = prog_vars,
                      "predictive" = pred_vars,
                      "both" = c(pred_vars, prog_vars))

  if (mnar) {
    covs <- c(covs, "U")
    data <- cbind(data, U)
  }

  # every combination of observed/missing over miss_vars, less the all-observed
  # and all-missing rows that ampute rejects
  if (length(miss_vars) > 1) {
    indicators <- expand.grid(rep(list(c(0, 1)), length(miss_vars)))
    indicators <- indicators[!apply(indicators, 1, function(x) all(x == 1)), ]
    colnames(indicators) <- miss_vars
    if (length(miss_vars) < length(covs)) {
      observed <- setdiff(covs, miss_vars)
      indicators[observed] <- 1
    }
    indicators <- indicators %>% select(all_of(covs))
    indicators <- indicators[!apply(indicators, 1, function(x) all(x == 0)), ]
  } else {
    indicators <- ifelse(covs == miss_vars, 1, 0)
    names(indicators) <- covs
  }

  # under MNAR only U drives the missingness, so it takes all the weight
  weights <- NULL
  if (mnar) {
    weights <- matrix(0, ncol = length(covs),
                      nrow = if (is.null(nrow(indicators))) 1 else nrow(indicators))
    weights[, length(covs)] <- 1
  }

  result <- ampute(data, prop = prop, patterns = indicators, weights = weights)

  data <- cbind(result$amp, keep)
  data %>% select(all_of(orig))
}

#' Apply a missing-data handling method
#'
#' @param data dataset containing NAs
#' @param method one of MISS_METHODS
#' @param n_imp number of imputations for method = "multiple_imputation"
#' @return list with `data`, plus `retained_indices` for the row-dropping methods
#'   and `ipw` for the IPW method. For multiple imputation `data` is a list of
#'   n_imp completed datasets.
handle_missingness <- function(data, method, n_imp = 50) {

  if (!any(is.na(data))) {
    message("No missing data found. Returning original dataset.")
    return(data)
  }
  method <- match.arg(method, MISS_METHODS)

  switch(method,
    "complete_cases" = {
      retained_indices <- complete.cases(data)
      complete_data <- data[retained_indices, ]
      message(paste("Complete case analysis: Removed", nrow(data) - nrow(complete_data)))
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
      # missForest wants binary columns as factors, so round-trip them
      df <- df %>%
        mutate(across(everything(), function(x) {
          unique_vals <- unique(x[!is.na(x)])
          bin <- length(unique_vals) == 2 && all(unique_vals %in% c(0, 1))
          if (bin) factor(x, levels = c(0, 1)) else x
        }))
      mf_imputed <- missForest(df)
      imputed_df <- as.data.frame(lapply(mf_imputed$ximp, function(x) {
        if (is.factor(x)) as.numeric(as.character(x)) else x
      }))
      message("Imputation with missForest complete")
      list(data = cbind(keep, imputed_df))
    },
    "regression" = {
      keep <- data %>% select(all_of(c("Y", "W")))
      df <- data %>% select(-all_of(c("Y", "W")))

      miss <- names(df)[sapply(df, function(x) any(is.na(x)))]
      complete <- setdiff(names(df), miss)
      imputed_df <- df
      for (var in miss) {
        fmla <- as.formula(paste(var, "~", paste(complete, collapse = " + ")))
        temp_imputed <- regressionImp(fmla, df)
        imputed_df[[var]] <- temp_imputed[[var]]
      }
      message("Imputation via regression complete")
      list(data = cbind(keep, imputed_df))
    },
    "missing_indicator" = {
      miss <- names(data)[sapply(data, function(x) any(is.na(x)))]
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

      miss <- names(df)[sapply(df, function(x) any(is.na(x)))]
      complete <- setdiff(names(df), miss)
      df$cc <- ifelse(complete.cases(df), 1, 0)

      fmla <- as.formula(paste("cc ~", paste(complete, collapse = " + ")))
      miss_lr <- glm(fmla, family = binomial, data = df)

      retained_indices <- complete.cases(df)
      complete_data <- data[retained_indices, ]
      ipw <- 1 / miss_lr$fitted.values[retained_indices]

      message(paste("IPW: removed", nrow(data) - nrow(complete_data), "observations"))
      message(paste("Final sample size:", nrow(complete_data)))
      list(data = complete_data, ipw = ipw, retained_indices = retained_indices)
    },
    "multiple_imputation" = {
      n_var <- ncol(data)
      predMat <- matrix(1, nrow = n_var, ncol = n_var)
      diag(predMat) <- 0
      predMat[, which(colnames(data) %in% c("Y", "W"))] <- 0

      imputation <- mice(data, m = n_imp, method = "rf", predictorMatrix = predMat)
      data_mi <- lapply(seq_len(n_imp), function(i) complete(imputation, i))
      message(paste0("Missing imputation completed - returning ", n_imp,
                     " imputed datasets"))
      list(data = data_mi)
    },
    "none" = {
      message("No missing data handling applied, original data set with missingness returned")
      list(data = data)
    })
}

#' Generate, amputate and handle in one call
#'
#' @param set scenario set, as for generate_scenario_data
#' @param n_imp number of imputations for method = "multiple_imputation"
generate_and_process_data <- function(scenario, n, set, return_truth = TRUE,
                                      type, prop, mech, method, n_imp = 50) {

  data_result <- generate_scenario_data(scenario, n, set,
                                        return_truth = return_truth, mech = mech)

  mech_norm <- sub("^AUX", "MNAR", mech)
  miss_dataset <- introduce_missingness(
    data_result$dataset, type, prop, mech,
    U = if (mech_norm %in% c("MNAR", "MNAR-Y")) data_result$truth$U else NULL)

  processed <- handle_missingness(miss_dataset, method, n_imp = n_imp)
  data_result$dataset <- processed$data
  data_result$missing_method <- method

  # the row-dropping methods must drop the same rows from the truth
  if (method %in% c("complete_cases", "IPW") && return_truth) {
    data_result$truth <- data_result$truth[processed$retained_indices, ]
  }
  if (method == "IPW") data_result$ipw <- processed$ipw

  data_result
}
