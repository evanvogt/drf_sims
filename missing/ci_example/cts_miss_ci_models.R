##########
# title: fit cate models and confidence intervals for mulitply-imputed data
##########
# Estimators in R/cate_models.R, bootstrap and MI pooling in R/bootstrap_ci.R.
#
# This study asks how to build a CATE confidence interval when the covariates
# were multiply imputed: fit every imputed dataset, then pool. The pooling
# offers three strategies side by side (pooled / mib / hybrid) - see
# combine_mi_ci().
#
# Profile "ci_mi" is the "ci" profile plus the nuisance row means, which this
# study kept and confidence_intervals/ did not.

source(here::here("R", "cate_models.R"))

run_all_cate_methods <- function(data, n_folds = 10, fmla_info = NULL,
                                 CI_boot = 200, CI_sf = 0.5, alpha = 0.05,
                                 num.threads = NULL, verbose_timing = FALSE) {
  cate_methods(data, n_folds = n_folds, sl_lib = NULL, fmla_info = fmla_info,
               family = gaussian(), profile = "ci_mi",
               ci = list(boot = CI_boot, sf = CI_sf, alpha = alpha),
               num.threads = num.threads, verbose_timing = verbose_timing)
}

#' Fit every imputed dataset and pool the results
#'
#' @param datalist list of completed datasets from multiple imputation
#' @param alpha two-sided level; now passed explicitly through to the pooling
#'   step, which previously read it from the global environment
#' @param num.threads grf thread count, forwarded to every regression_forest()/
#'   causal_forest() call inside cate_methods(). NULL preserves the prior
#'   behaviour (grf's own default). Added for cts_miss_ci_profile.R to be able
#'   to control it; cts_miss_ci_analysis.R now also forwards it from its own
#'   trailing commandArgs.
mi_boot <- function(datalist, n_folds = 10, fmla_info = NULL,
                    CI_boot = 200, CI_sf = 0.5, alpha = 0.05,
                    num.threads = NULL) {

  results_list <- future_map(datalist, function(data) {
    run_all_cate_methods(
      data = data,
      n_folds = n_folds,
      fmla_info = fmla_info,
      CI_boot = CI_boot,
      CI_sf = CI_sf,
      alpha = alpha,
      num.threads = num.threads
    )
  }, .options = furrr_options(seed = TRUE))

  results <- list()
  results$causal_forest <- combine_mi_ci(results_list, "causal_forest", alpha)
  results$dr_random_forest <- combine_mi_ci(results_list, "dr_random_forest", alpha)
  if (!is.null(fmla_info)) {
    results$dr_oracle <- combine_mi_ci(results_list, "dr_oracle", alpha)
  }
  results$dr_semi_oracle <- combine_mi_ci(results_list, "dr_semi_oracle", alpha)

  results
}
