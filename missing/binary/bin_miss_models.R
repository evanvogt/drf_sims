##########
# title: CATE model fitting functions - missing covariates and binary outcome
##########
# As missing/continuous/cts_miss_models.R, with a binomial outcome family.
# Estimators in R/cate_models.R.
#
# The oracle formula from R/dgm_scenarios.R returns the risk itself, so the
# oracle arm needs no link.
#
# Bug M lived here: this file passed oracle_link = "identity" after the binary
# missing-data table lost its plogis(...)-wrapped oracle strings (6b06db3), so
# dr_oracle was handed log-odds as outcome predictions for a 0/1 outcome. No
# finished result was affected - see missing/binary/README.md. Since the
# risk-difference DGM there is no oracle_link argument left to get wrong.

source(here::here("R", "cate_models.R"))

run_all_cate_methods <- function(data, n_folds = 10, sl_lib = NULL,
                                 fmla_info = NULL, ipw = NULL,
                                 num.threads = NULL, verbose_timing = FALSE) {
  cate_methods(data, n_folds = n_folds, sl_lib = sl_lib, fmla_info = fmla_info,
               family = binomial(),
               ipw = ipw, profile = "missing",
               num.threads = num.threads, verbose_timing = verbose_timing)
}
