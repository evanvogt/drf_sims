##########
# title: CATE model fitting functions - missing covariates and binary outcome
##########
# As missing/continuous/cts_miss_models.R, with a binomial outcome family.
# Estimators in R/cate_models.R.
#
# oracle_link = "logit", as binary/bin_models.R: the oracle formula from
# R/dgm_scenarios.R is a LINEAR PREDICTOR, so the model applies plogis.
#
# This was "identity" until bug M. The binary missing-data table used to carry
# its own plogis(...)-wrapped oracle strings; deleting that legacy table (6b06db3)
# rebuilt binary_missing_fixed on continuous_missing's plain strings and left this
# link alone, so dr_oracle was handed log-odds as outcome predictions for a 0/1
# outcome. No finished result was affected - see missing/binary/README.md.

source(here::here("R", "cate_models.R"))

run_all_cate_methods <- function(data, n_folds = 10, sl_lib = NULL,
                                 fmla_info = NULL, ipw = NULL,
                                 num.threads = NULL, verbose_timing = FALSE) {
  cate_methods(data, n_folds = n_folds, sl_lib = sl_lib, fmla_info = fmla_info,
               family = binomial(), oracle_link = "logit",
               ipw = ipw, profile = "missing",
               num.threads = num.threads, verbose_timing = verbose_timing)
}
