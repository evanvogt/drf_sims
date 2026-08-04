###############
# title: data generating process - binary outcomes, all scenarios
###############
# Scenario tables and generator in R/dgm_scenarios.R. The same ten scenarios as
# the continuous study on a logit scale, differing in the coefficient table and
# in scenario 10's treatment effect (exp(X4) rather than exp(-abs(X4))).
#
# The oracle formula here is a LINEAR PREDICTOR and bin_models.R applies plogis
# (oracle_link = "logit"). missing/binary/ uses the opposite convention.

source(here::here("R", "dgm_scenarios.R"))

generate_binary_scenario_data <- function(scenario, n, return_truth = TRUE,
                                          seed = NULL) {
  generate_scenario_data(scenario, n, set = "binary",
                         return_truth = return_truth, seed = seed)
}

get_binary_oracle_info <- function(scenario, bW) {
  get_oracle_info(scenario, bW, set = "binary")
}
