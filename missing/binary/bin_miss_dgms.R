###############
# title: binary DGMs with missingness and imputations
###############
# Scenarios in R/dgm_scenarios.R, missingness machinery in R/missingness.R.
#
# Previously carried the continuous coefficient table, wrong t-test calibration,
# and un-plogis'd truth (three related defects, fixed together).

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
