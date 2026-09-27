###############
# title: data generating process - binary outcomes, all scenarios
###############
# Scenario tables and generator in R/dgm_scenarios.R. The same ten scenarios as
# the continuous study, with the treatment effect on the RISK-DIFFERENCE scale:
# P(Y = 1) = m0(x) + W * tau(x), the control risk m0 bounded in [0.34, 0.70],
# X4/X5 entering tau through tanh, and the modifiers sign-reversed and rescaled
# (RD_SCALE). See README.md, "Outcome model and bW calibration".
#
# The oracle formula here returns the risk itself.

source(here::here("R", "dgm_scenarios.R"))

generate_binary_scenario_data <- function(scenario, n, return_truth = TRUE,
                                          seed = NULL) {
  generate_scenario_data(scenario, n, set = "binary",
                         return_truth = return_truth, seed = seed)
}

get_binary_oracle_info <- function(scenario, bW) {
  get_oracle_info(scenario, bW, set = "binary")
}
