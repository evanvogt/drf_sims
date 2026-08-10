###############
# title: data generating process - interim-analysis validation study
###############
# The scenario tables and the generator live in R/dgm_scenarios.R. This file
# names this study's slice of them - the same continuous-outcome slice that
# continuous/cts_dgms.R names, since this study fits on chunks of the same
# scenario. It replaces the standalone fork that used to live in
# cts_dgm_validation.R.

source(here::here("R", "dgm_scenarios.R"))

generate_continuous_scenario_data <- function(scenario, n, return_truth = TRUE,
                                              seed = NULL) {
  generate_scenario_data(scenario, n, set = "continuous",
                         return_truth = return_truth, seed = seed)
}

get_continuous_oracle_info <- function(scenario, bW) {
  get_oracle_info(scenario, bW, set = "continuous")
}
