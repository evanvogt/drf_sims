##########
# title: CATE estimation functions - binary outcome
##########
# The estimators live in R/cate_models.R. This file is the binary study's
# configuration of them, and the outcome type is the entire difference from the
# continuous study across what used to be 438 near-identical lines:
#   family = binomial()      SuperLearner outcome model + method.NNloglik
#   oracle_link = "logit"    bin_dgms.R's get_binary_oracle_info returns a linear
#                            predictor, so plogis is applied by the model code.
#                            (missing/binary/ uses the opposite convention.)

source(here::here("R", "cate_models.R"))

run_all_cate_methods <- function(data, n_folds = 10, sl_lib = NULL, fmla_info = NULL,
                                 num.threads = NULL, verbose_timing = FALSE) {
  cate_methods(data, n_folds = n_folds, sl_lib = sl_lib, fmla_info = fmla_info,
               family = binomial(), oracle_link = "logit", profile = "base",
               num.threads = num.threads, verbose_timing = verbose_timing)
}
