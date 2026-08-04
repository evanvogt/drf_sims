##########
# title: CATE model fitting functions - missing covariates and binary outcome
##########
# As missing/continuous/cts_miss_models.R, with a binomial outcome family.
# Estimators in R/cate_models.R.
#
# oracle_link = "identity" - NOT a mistake, and NOT the same as binary/. This
# study's bin_miss_dgms.R bakes plogis(...) into the oracle formula string
# itself, so the model must not apply the link a second time. binary/bin_dgms.R
# uses the opposite convention. The two give identical oracle values; only the
# division of labour differs. See the header of R/cate_models.R.

source(here::here("R", "cate_models.R"))

run_all_cate_methods <- function(data, n_folds = 10, sl_lib = NULL,
                                 fmla_info = NULL, ipw = NULL) {
  cate_methods(data, n_folds = n_folds, sl_lib = sl_lib, fmla_info = fmla_info,
               family = binomial(), oracle_link = "identity",
               ipw = ipw, profile = "missing")
}
