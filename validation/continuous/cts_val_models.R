##########
# title: CATE estimation for the interim-analysis validation study - continuous
##########
# The estimators, both importance measures, the interaction tests and the four
# chunk comparisons are shared with binary/ and live in validation/val_common.R.
# This file is the continuous arm's configuration of them, and the outcome type
# is the entire difference from binary/bin_val_models.R:
#   family = gaussian()      the DR SuperLearner's outcome model

source(here::here("validation", "val_common.R"))

#' Fit the three estimators and both importance measures on one trial chunk
#'
#' fit_val_methods() (validation/val_common.R) with a continuous outcome - see
#' there for the arguments.
run_all_cate_methods <- function(data, n_folds = 10, num.threads = NULL,
                                 verbose_timing = FALSE,
                                 sl_lib = sl_libraries(nrow(data))) {
  fit_val_methods(data, n_folds = n_folds, num.threads = num.threads,
                  verbose_timing = verbose_timing, sl_lib = sl_lib,
                  family = gaussian())
}
