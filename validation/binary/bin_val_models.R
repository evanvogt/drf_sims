##########
# title: CATE estimation for the interim-analysis validation study - binary
##########
# The estimators, both importance measures, the interaction tests and the four
# chunk comparisons are shared with continuous/ and live in
# validation/val_common.R. This file is the binary arm's configuration of them,
# and the outcome type is the entire difference from
# continuous/cts_val_models.R:
#   family = binomial()      the DR SuperLearner's outcome model
#                            (method.NNloglik)
# The interaction tests' robust = TRUE is the analysis script's to pass, not
# this file's - see bin_val_analysis.R.

source(here::here("validation", "val_common.R"))

#' Fit the three estimators and both importance measures on one trial chunk
#'
#' fit_val_methods() (validation/val_common.R) with a binary outcome - see there
#' for the arguments.
run_all_cate_methods <- function(data, n_folds = 10, num.threads = NULL,
                                 verbose_timing = FALSE,
                                 sl_lib = sl_libraries(nrow(data))) {
  fit_val_methods(data, n_folds = n_folds, num.threads = num.threads,
                  verbose_timing = verbose_timing, sl_lib = sl_lib,
                  family = binomial())
}
