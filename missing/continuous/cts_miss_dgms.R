###############
# title: continuous DGMs with missingness and imputations
###############
# Scenarios in R/dgm_scenarios.R, missingness machinery in R/missingness.R.
#
# Five renumbered scenarios, NOT a subset of the main study's ten: scenario k
# here corresponds to 1, 2, 4, 8, 9 there.

source(here::here("R", "missingness.R"))

generate_continuous_scenario_data <- function(scenario, n, mech, return_truth = TRUE) {
  generate_scenario_data(scenario, n, set = "continuous_missing",
                         return_truth = return_truth, mech = mech)
}

introduce_missingness_continuous <- function(data, type, prop, mech, U = NULL) {
  introduce_missingness(data, type, prop, mech, U = U)
}

handle_missingness_continuous <- function(data, method) {
  handle_missingness(data, method, n_imp = 50)
}

generate_and_process_continuous_data <- function(scenario, n, return_truth = TRUE,
                                                 type, prop, mech, method) {
  generate_and_process_data(scenario, n, set = "continuous_missing",
                            return_truth = return_truth, type = type,
                            prop = prop, mech = mech, method = method,
                            n_imp = 50)
}

get_continuous_oracle_info <- function(scenario, bW) {
  get_oracle_info(scenario, bW, set = "continuous_missing")
}
