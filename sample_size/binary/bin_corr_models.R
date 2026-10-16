##########
# title: CATE estimation functions - correlated binary outcome
##########
# The estimators live in R/cate_models.R. This file is the correlated binary
# study's configuration of them, and the outcome type is the entire difference
# from the continuous one (cts_corr_models.R):
#   family = binomial()      SuperLearner outcome model + method.NNloglik
# The oracle formula from the binary DGM returns the risk itself, so the oracle
# arm needs no link. It was sample_size/binary/bin_models.R until that study was
# retired (2026-10-06); the function is unchanged.

source(here::here("R", "cate_models.R"))

run_all_cate_methods <- function(data, n_folds = 10, sl_lib = NULL, fmla_info = NULL,
                                 num.threads = NULL, verbose_timing = FALSE) {
  cate_methods(data, n_folds = n_folds, sl_lib = sl_lib, fmla_info = fmla_info,
               family = binomial(), profile = "base",
               num.threads = num.threads, verbose_timing = verbose_timing)
}
