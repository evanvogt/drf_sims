##########
# title: CATE model fitting functions - missing covariates and continuous outcome
##########
# The estimators live in R/cate_models.R. This file is the missing-data
# continuous study's configuration of them: the "missing" orchestration profile,
# which skips the arms needing a complete covariate matrix when X still has NAs,
# and threads `ipw` through as grf sample.weights / SuperLearner obsWeights.
#
# combine_mi (Rubin's rules across imputations) also lives in R/cate_models.R.

source(here::here("R", "cate_models.R"))

run_all_cate_methods <- function(data, n_folds = 10, sl_lib = NULL,
                                 fmla_info = NULL, ipw = NULL) {
  cate_methods(data, n_folds = n_folds, sl_lib = sl_lib, fmla_info = fmla_info,
               family = gaussian(), ipw = ipw, profile = "missing")
}
