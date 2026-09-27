##########
# title: CATE model fitting with confidence intervals - cts outcome
##########
# The estimators live in R/cate_models.R and the bootstrap in R/bootstrap_ci.R.
# This file is the continuous CI study's configuration of them: the "ci"
# orchestration profile, which adds half-sample bootstrap intervals and
# causal-forest variance estimates, and drops the SuperLearner arm and the
# post-estimation BLP / independence tests.

source(here::here("R", "cate_models.R"))

run_all_cate_methods <- function(data, n_folds = 10, fmla_info = NULL,
                                 CI_boot = 200, CI_sf = 0.5, alpha = 0.05,
                                 Z_query = NULL) {
  cate_methods(data, n_folds = n_folds, sl_lib = NULL, fmla_info = fmla_info,
               family = gaussian(), profile = "ci",
               ci = list(boot = CI_boot, sf = CI_sf, alpha = alpha),
               Z_query = Z_query)
}
