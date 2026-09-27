##########
# title: CATE model fitting with confidence intervals - binary outcome
##########
# The estimators live in R/cate_models.R and the bootstrap in R/bootstrap_ci.R.
# As the continuous CI study, but with a binomial outcome family. The oracle
# formula from bin_ci_dgms.R returns the risk itself, so the oracle arm needs no
# link. (Bug A - the DGM carrying the continuous coefficient table - is fixed in
# the DGM; the models were never affected.)

source(here::here("R", "cate_models.R"))

run_all_cate_methods <- function(data, n_folds = 10, fmla_info = NULL,
                                 CI_boot = 200, CI_sf = 0.5, alpha = 0.05,
                                 Z_query = NULL) {
  cate_methods(data, n_folds = n_folds, sl_lib = NULL, fmla_info = fmla_info,
               family = binomial(), profile = "ci",
               ci = list(boot = CI_boot, sf = CI_sf, alpha = alpha),
               Z_query = Z_query)
}
