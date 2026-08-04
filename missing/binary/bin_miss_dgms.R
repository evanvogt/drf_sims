###############
# title: binary DGMs with missingness and imputations
###############
# Scenarios in R/dgm_scenarios.R, missingness machinery in R/missingness.R.
#
# TWO KNOWN BUGS, both gated by flags in R/dgm_scenarios.R:
#
#  LEGACY_BIN_MISS_PARAMS - this study carries the CONTINUOUS coefficient table
#    (b0 = c(0.4, 0.2, 0.4, 1, 0.4), b1 = -0.05, b2 = c(2, 2, 2, 2, 1)) rather
#    than the binary one, the same copy-paste as bug A in
#    confidence_intervals/binary.
#
#  LEGACY_BIN_MISS_TRUTH - truth was computed as p0 = b0 + b1*X1 + b2*X2 with no
#    plogis, so truth$tau is a difference in LOG-ODDS, while the outcome is
#    rbinom(n, 1, plogis(lp)) and every estimator targets a RISK DIFFERENCE.
#    bin_miss_metrics.R compares the two directly.
#
# Flipping either invalidates everything under ../results/missing/binary.
#
# The oracle formula here already contains plogis(...), so bin_miss_models.R must
# pass oracle_link = "identity" - the opposite convention to binary/bin_dgms.R.

source(here::here("R", "missingness.R"))

generate_binary_scenario_data <- function(scenario, n, mech, return_truth = TRUE) {
  generate_scenario_data(scenario, n, set = "binary_missing",
                         return_truth = return_truth, mech = mech)
}

introduce_missingness_binary <- function(data, type, prop, mech, U = NULL) {
  introduce_missingness(data, type, prop, mech, U = U)
}

handle_missingness_binary <- function(data, method) {
  handle_missingness(data, method, n_imp = 50)
}

generate_and_process_binary_data <- function(scenario, n, return_truth = TRUE,
                                             type, prop, mech, method) {
  generate_and_process_data(scenario, n, set = "binary_missing",
                            return_truth = return_truth, type = type,
                            prop = prop, mech = mech, method = method,
                            n_imp = 50)
}

get_binary_oracle_info <- function(scenario, bW) {
  get_oracle_info(scenario, bW, set = "binary_missing")
}
