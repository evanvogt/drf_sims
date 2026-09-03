###############
# title: data generating mechanisms for CIs - bin outcome
###############
# Previously ran on the continuous coefficient table (bug A); fixed — now always uses the binary table.

source(here::here("R", "dgm_scenarios.R"))

generate_binary_scenario_data <- function(scenario, n, return_truth = TRUE,
                                          seed = NULL) {
  generate_scenario_data(scenario, n, set = "binary_ci",
                         return_truth = return_truth, seed = seed)
}

get_binary_oracle_info <- function(scenario, bW) {
  get_oracle_info(scenario, bW, set = "binary_ci")
}
