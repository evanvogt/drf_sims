###############
# title: continuous DGMs with missingness - CI example
###############
# Same scenarios and machinery as missing/continuous/, reached through
# R/dgm_scenarios.R and R/missingness.R.
#
# Its grid runs MAR only, so the unobserved-U branches never fire here. It uses
# continuous_missing's correlated covariates (X01-X03 auxiliary).

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
