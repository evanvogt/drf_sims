##########
# title: shared missingness generation and handling
##########
# One copy of what was duplicated across missing/continuous/, missing/binary/ and
# missing/ci_example/. The three differed only in the function-name suffix
# (_continuous / _binary), the MNAR-vs-AUX spelling, and the number of multiple
# imputations - none of which is about the outcome type, so the suffixes are gone
# and the imputation count is an argument. The mechanisms are MISS_MECHS
# (R/dgm_scenarios.R) since 2026-09-28; the old names are rejected.

require(dplyr)
require(mice)
require(missForest)
require(VIM)

source(here::here("R", "dgm_scenarios.R"))

# columns that are never amputated: the outcome, the treatment, and X01-X05 -
# outside the outcome model, and in the missing-data sets X01-X03 are the
# auxiliaries of the amputed covariates (R/dgm_scenarios.R)
NEVER_MISSING <- c("Y", "W", "X01", "X02", "X03", "X04", "X05")
PROGNOSTIC_VARS <- c("X1", "X2")

MISS_METHODS <- c("complete_cases", "mean_imputation", "missforest", "regression",
                  "missing_indicator", "IPW", "multiple_imputation", "none")

#' Which covariates an amputation left missing
#'
#' Saved with every missing-data run (since 2026-09-29) so the metrics can
#' score each method on the complete units, the incomplete units, and the
#' complete-unit subset every method shares (R/metrics.R,
#' cate_metrics_split()). The truth is tau at the unamputed covariates, which no
#' method can recover for an incomplete unit, so without the split the
#' imputation arms carry an error floor the row-dropping arms never face.
#'
#' @param amputed introduce_missingness()'s output
#' @return n x k logical matrix over the amputable covariates (X1-X5)
amputation_mask <- function(amputed) {
  is.na(as.matrix(amputed[, setdiff(names(amputed), NEVER_MISSING), drop = FALSE]))
}

#' Apply an imputation separately within each treatment arm
#'
#' Imputing within arm keeps every W x X interaction in the imputation model,
#' which a pooled model with main effects would flatten - the point of
#' imputing for a CATE. W is constant within an arm, so `impute` is handed the
#' arm's rows without it and must return them (or a list of completed copies)
#' with the same columns; W is put back and the rows returned in the original
#' order, since cate_methods() reads covariates by position. An arm with
#' nothing missing is returned as it is.
#'
#' Every arm present in W is imputed on its own, in increasing order of W, so
#' a multi-arm W (e.g. 0, 1, 2) works as well; for a binary W that is the 0, 1
#' order this always had.
#'
#' @param data dataset with Y, W and covariates, containing NAs
#' @param impute function(df) returning a completed data.frame, or a list of
#'   them (multiple imputation)
#' @param n_out 1, or the number of completed datasets impute returns
#' @return a data.frame, or a list of n_out data.frames
impute_by_arm <- function(data, impute, n_out = 1) {
  out <- rep(list(data), n_out)
  for (arm in sort(unique(data$W))) {
    rows <- which(data$W == arm)
    df <- data[rows, setdiff(names(data), "W"), drop = FALSE]
    if (!anyNA(df)) next
    done <- impute(df)
    if (is.data.frame(done)) done <- list(done)
    for (i in seq_len(n_out)) {
      cols <- names(done[[i]])
      out[[i]][rows, cols] <- done[[i]][, cols]
    }
  }
  if (n_out == 1) out[[1]] else out
}

#' Introduce missingness into a simulated dataset
#'
#' Builds an mice::ampute pattern matrix over the covariates being amputated and
#' applies it. Under MAR the missingness depends on the observed covariates;
#' under MNAR-Y0 / MNAR-tau it is driven entirely by the unobserved U, which is
#' what the weight matrix below encodes. The two MNAR mechanisms amputate
#' identically; they differ only in where the generator puts U in the outcome
#' (R/dgm_scenarios.R).
#'
#' @param data simulated dataset
#' @param type which covariates to amputate: "prognostic", "predictive" or "both"
#' @param prop proportion of missingness, in (0, 1)
#' @param mech one of MISS_MECHS: "MAR", "MNAR-Y0" or "MNAR-tau"
#' @param U the unobserved variable, required for the MNAR mechanisms
introduce_missingness <- function(data, type, prop, mech, U = NULL) {

  check_mech(mech)

  if (!type %in% c("prognostic", "predictive", "both")) {
    stop("type must be 'prognostic', 'predictive', or 'both'")
  }
  if (prop < 0 || prop > 1) stop("miss_prop must be between 0 and 1")
  mnar <- mech %in% MNAR_MECHS
  if (mnar && is.null(U)) {
    stop("unobserved variable U required for MNAR missingness generation")
  }

  orig <- colnames(data)
  keep <- data %>% select(all_of(NEVER_MISSING))
  data <- data %>% select(-all_of(NEVER_MISSING))
  covs <- colnames(data)

  # every scenario draws X1-X5 (R/dgm_scenarios.R), so "predictive" is X3-X5
  # whether or not the scenario's treatment effect uses them, and "both" - what
  # every missing-data grid runs - amputates all of X1-X5
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

  # mice::ampute (3.19.0, ampute.continuous()) has an edge case: a pattern with
  # exactly one candidate whose weighted sum score is 0 gets the response
  # indicator R <- 0, a scalar, which recycles over EVERY row - so that
  # pattern's variables go missing for the whole dataset, and rows turn up with
  # missingness patterns that were never specified (all-missing, unions of
  # patterns). Found 2026-09-29: 6 of 40 MAR draws at n = 120 (about 4
  # candidates per pattern). At the studies' n = 500 a pattern has about 17
  # candidates and a single-candidate pattern is rare (none in 3000 MAR and
  # 1000 MNAR-Y0 draws of scenario 2), but
  # guard anyway: redraw until every incomplete row carries one of the
  # specified patterns. A valid first draw - every one at n = 500 so far - is
  # returned unchanged, with the same RNG consumption as before the guard.
  pattern_key <- function(m) apply(m, 1, paste, collapse = "")
  allowed <- pattern_key(matrix(as.matrix(indicators) == 0, ncol = length(covs)))
  for (attempt in 1:20) {
    result <- ampute(data, prop = prop, patterns = indicators, weights = weights)
    na <- is.na(as.matrix(result$amp[, covs, drop = FALSE]))
    incomplete <- rowSums(na) > 0
    if (all(pattern_key(na[incomplete, , drop = FALSE]) %in% allowed)) break
    if (attempt == 20) {
      stop("mice::ampute produced unspecified missingness patterns 20 times running",
           call. = FALSE)
    }
    warning("mice::ampute produced unspecified missingness patterns (a pattern ",
            "with a single candidate) - redrawing the amputation", call. = FALSE)
  }

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
      # within each arm, with Y as a predictor (since 2026-09-29; before then Y
      # and W were both left out, which flattened the HTE - see impute_by_arm)
      imputed <- impute_by_arm(data, function(df) {
        # missForest wants binary columns as factors (Y too, for a binary
        # outcome), so round-trip them
        df <- df %>%
          mutate(across(everything(), function(x) {
            unique_vals <- unique(x[!is.na(x)])
            bin <- length(unique_vals) == 2 && all(unique_vals %in% c(0, 1))
            if (bin) factor(x, levels = c(0, 1)) else x
          }))
        mf_imputed <- missForest(df)
        as.data.frame(lapply(mf_imputed$ximp, function(x) {
          if (is.factor(x)) as.numeric(as.character(x)) else x
        }))
      })
      message("Imputation with missForest complete")
      list(data = imputed)
    },
    "regression" = {
      # single deterministic regression imputation, within each arm, each
      # incomplete covariate regressed (lm) on the fully observed columns: the
      # always-observed covariates and Y (since 2026-09-30; before then one
      # pooled model on the covariates alone, Y and W excluded, which flattened
      # the HTE - see impute_by_arm). The predictor set is fixed from the whole
      # dataset, so both arms fit the same model even if an arm happens to
      # observe every value of some covariate.
      miss <- names(data)[sapply(data, anyNA)]
      complete <- setdiff(names(data), c(miss, "W"))
      imputed <- impute_by_arm(data, function(df) {
        for (var in miss) {
          fmla <- as.formula(paste(var, "~", paste(complete, collapse = " + ")))
          df[[var]] <- regressionImp(fmla, df)[[var]]
        }
        df
      })
      message("Imputation via regression complete")
      list(data = imputed)
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
      # P(complete) by logistic regression on the fully observed covariates
      # plus W * Y, fit on every unit (since 2026-09-30; before then Y and W
      # were excluded, so under MNAR the weights were near constant and IPW
      # reproduced complete cases). W alone cannot predict the missingness of
      # a baseline covariate; W * Y lets its relation to Y differ by arm.
      covs <- setdiff(names(data), c("Y", "W"))
      miss <- covs[sapply(data[covs], anyNA)]
      complete <- setdiff(covs, miss)
      retained_indices <- complete.cases(data[covs])
      df <- data
      df$cc <- as.numeric(retained_indices)

      fmla <- as.formula(paste("cc ~", paste(c(complete, "W * Y"), collapse = " + ")))
      miss_lr <- glm(fmla, family = binomial, data = df)

      complete_data <- data[retained_indices, ]
      # stabilised by the marginal P(complete), so the weights average about 1
      # over the complete cases. A constant rescaling, but not a no-op: grf
      # (2.5.0) forests change with the weights' scale - constant weights of 1
      # reproduce the unweighted forest exactly, constant weights of 1.43
      # (= 1 / 0.7, the unstabilised scale) do not - so stabilising keeps the
      # weighted forests on the same footing as every other arm's
      # (checked 2026-09-30)
      ipw <- mean(retained_indices) / miss_lr$fitted.values[retained_indices]

      message(paste("IPW: removed", nrow(data) - nrow(complete_data), "observations"))
      message(paste("Final sample size:", nrow(complete_data)))
      list(data = complete_data, ipw = ipw, retained_indices = retained_indices)
    },
    "multiple_imputation" = {
      # within each arm, every other column - Y included - predicting each
      # incomplete one (since 2026-09-29; before then one pooled model with Y
      # and W both excluded, which flattened the HTE - see impute_by_arm)
      data_mi <- impute_by_arm(data, function(df) {
        n_var <- ncol(df)
        predMat <- matrix(1, nrow = n_var, ncol = n_var,
                          dimnames = list(names(df), names(df)))
        diag(predMat) <- 0
        imputation <- mice(df, m = n_imp, method = "rf", predictorMatrix = predMat)
        lapply(seq_len(n_imp), function(i) complete(imputation, i))
      }, n_out = n_imp)
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

  miss_dataset <- introduce_missingness(
    data_result$dataset, type, prop, mech,
    U = if (mech %in% MNAR_MECHS) data_result$truth$U else NULL)
  data_result$miss_mask <- amputation_mask(miss_dataset)

  processed <- handle_missingness(miss_dataset, method, n_imp = n_imp)
  data_result$dataset <- processed$data
  data_result$missing_method <- method

  # which of the n units the analysed dataset keeps: the complete cases for
  # the row-dropping methods, every unit otherwise
  data_result$retained_indices <- if (method %in% c("complete_cases", "IPW")) {
    processed$retained_indices
  } else {
    rep(TRUE, n)
  }

  # the row-dropping methods must drop the same rows from the truth
  if (method %in% c("complete_cases", "IPW") && return_truth) {
    data_result$truth <- data_result$truth[processed$retained_indices, ]
  }
  if (method == "IPW") data_result$ipw <- processed$ipw

  data_result
}

#' E[U_term | complete] and E[U_term | incomplete] under the MNAR mechanisms
#'
#' Under MNAR-Y0 / MNAR-tau the amputation selects on U alone, so the units
#' left complete have lower U than average and the incomplete ones higher.
#' Under MNAR-tau U's term is in the treatment effect, so the CATE given the
#' covariates AND whether a unit is complete is tau(X) plus one of these two
#' constants - the secondary truth cate_metrics_split() scores against
#' (missing/ADEMP.md, "Estimands"). They are population constants of the DGM,
#' estimated once here by a large simulation; the amputation depends on U
#' alone, so any scenario gives them. The caller's RNG state is restored.
#'
#' @param set a missing-data scenario set (one with a u_expr)
#' @param type,prop the amputation, as introduce_missingness()
#' @param n simulated sample size
#' @return c(complete = , incomplete = )
u_term_by_completeness <- function(set, type = "both", prop = 0.3, n = 1e5,
                                   seed = 20260929) {
  if (exists(".Random.seed", envir = globalenv())) {
    old_seed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()))
  }
  set.seed(seed)
  params <- resolve_set(set)
  # bW calibrated at the studies' n = 500, as they run (it does not move U or
  # the amputation, but keeps every binary risk where the studies' are)
  g <- generate_scenario_data(2, n, set, mech = "MNAR-tau", calib_n = 500)
  amputed <- suppressWarnings(
    introduce_missingness(g$dataset, type, prop, "MNAR-tau", U = g$truth$U))
  complete <- rowSums(amputation_mask(amputed)) == 0
  u_term <- eval(parse(text = params$u_expr[1]),
                 envir = list(bU = params$bU[1], U = g$truth$U))
  c(complete = mean(u_term[complete]), incomplete = mean(u_term[!complete]))
}

#' Generate the complete-data reference arm, with the mask it would have had
#'
#' The complete_data arm analyses the unamputed data, but its metrics need the
#' same complete / incomplete split as the other arms (amputation_mask()). The
#' seed stream is set by run alone, so generating and then amputating here
#' consumes exactly the draws generate_and_process_data() does, and the mask
#' is the one every other method of this (scenario, mechanism, run) sees. The
#' amputed copy is discarded; only its mask is kept.
#'
#' @inheritParams generate_and_process_data
generate_reference_data <- function(scenario, n, set, type, prop, mech) {
  data_result <- generate_scenario_data(scenario, n, set, return_truth = TRUE,
                                        mech = mech)
  amputed <- introduce_missingness(
    data_result$dataset, type, prop, mech,
    U = if (mech %in% MNAR_MECHS) data_result$truth$U else NULL)
  data_result$miss_mask <- amputation_mask(amputed)
  data_result$retained_indices <- rep(TRUE, n)
  data_result$missing_method <- "complete_data"
  data_result
}
