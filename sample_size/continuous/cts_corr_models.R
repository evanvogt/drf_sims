##########
# title: CATE estimation functions - correlated continuous outcome
##########
# The estimators live in R/cate_models.R. This file is the correlated continuous
# study's configuration of them: a Gaussian outcome on the "base" orchestration
# profile (post-estimation BLP and independence tests on, no bootstrap, no IPW
# weights). It was sample_size/continuous/cts_models.R until that study was
# retired (2026-10-06); the function is unchanged.

source(here::here("R", "cate_models.R"))

run_all_cate_methods <- function(data, n_folds = 10, sl_lib = NULL, fmla_info = NULL,
                                 num.threads = NULL, verbose_timing = FALSE) {
  cate_methods(data, n_folds = n_folds, sl_lib = sl_lib, fmla_info = fmla_info,
               family = gaussian(), profile = "base",
               num.threads = num.threads, verbose_timing = verbose_timing)
}
