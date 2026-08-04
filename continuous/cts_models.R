##########
# title: CATE estimation functions - continuous outcome
##########
# The estimators live in R/cate_models.R. This file is the continuous study's
# configuration of them: a Gaussian outcome on the "base" orchestration profile
# (post-estimation BLP and independence tests on, no bootstrap, no IPW weights).

source(here::here("R", "cate_models.R"))

run_all_cate_methods <- function(data, n_folds = 10, sl_lib = NULL, fmla_info = NULL) {
  cate_methods(data, n_folds = n_folds, sl_lib = sl_lib, fmla_info = fmla_info,
               family = gaussian(), profile = "base")
}
