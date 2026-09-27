###############
# title: data generating mechanisms for CIs - cts outcome
###############
# This file used to carry the note "just a straight copy of the cts dgm for
# now", and it was - byte for byte. It now names the same scenario set as
# continuous/cts_dgms.R rather than duplicating it.

source(here::here("R", "dgm_scenarios.R"))

generate_continuous_scenario_data <- function(scenario, n, return_truth = TRUE,
                                              seed = NULL) {
  generate_scenario_data(scenario, n, set = "continuous",
                         return_truth = return_truth, seed = seed)
}

get_continuous_oracle_info <- function(scenario, bW) {
  get_oracle_info(scenario, bW, set = "continuous")
}
