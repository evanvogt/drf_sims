##########
# title: CATE model fitting with confidence intervals - correlated binary outcome
##########
# The estimators live in R/cate_models.R and the bootstrap in R/bootstrap_ci.R.
# As cts_corr_ci_models.R, but with a binomial outcome family; used by the
# correlated binary CI and optimal_sf studies. The oracle formula from the
# binary DGM returns the risk itself, so the oracle arm needs no link. It was
# sample_size/confidence_intervals/binary/bin_ci_models.R until that study was
# retired (2026-10-06); the function is unchanged.

source(here::here("R", "cate_models.R"))

run_all_cate_methods <- function(data, n_folds = 10, fmla_info = NULL,
                                 CI_boot = 200, CI_sf = 0.5, alpha = 0.05,
                                 Z_query = NULL) {
  cate_methods(data, n_folds = n_folds, sl_lib = NULL, fmla_info = fmla_info,
               family = binomial(), profile = "ci",
               ci = list(boot = CI_boot, sf = CI_sf, alpha = alpha),
               Z_query = Z_query)
}
