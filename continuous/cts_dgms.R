###############
# title: data generating process - continuous outcomes, all scenarios
###############
# The scenario tables and the generator live in R/dgm_scenarios.R. This file
# names the continuous study's slice of them.

source(here::here("R", "dgm_scenarios.R"))

generate_continuous_scenario_data <- function(scenario, n, return_truth = TRUE,
                                              seed = NULL) {
  generate_scenario_data(scenario, n, set = "continuous",
                         return_truth = return_truth, seed = seed)
}

get_continuous_oracle_info <- function(scenario, bW) {
  get_oracle_info(scenario, bW, set = "continuous")
}
