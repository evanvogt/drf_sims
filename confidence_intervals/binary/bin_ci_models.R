##########
# title: CATE model fitting with confidence intervals - binary outcome
##########
# The estimators live in R/cate_models.R and the bootstrap in R/bootstrap_ci.R.
# As the continuous CI study, but with a binomial outcome family.
#
# NOTE: this study's DGM (bin_ci_dgms.R) currently carries the CONTINUOUS
# coefficient table - see bug A in the repo README. The models are unaffected;
# the DGM is fixed separately.

source(here::here("R", "cate_models.R"))

run_all_cate_methods <- function(data, n_folds = 10, fmla_info = NULL,
                                 CI_boot = 200, CI_sf = 0.5, alpha = 0.05) {
  cate_methods(data, n_folds = n_folds, sl_lib = NULL, fmla_info = fmla_info,
               family = binomial(), oracle_link = "logit", profile = "ci",
               ci = list(boot = CI_boot, sf = CI_sf, alpha = alpha))
}
