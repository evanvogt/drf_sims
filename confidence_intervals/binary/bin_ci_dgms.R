###############
# title: data generating mechanisms for CIs - bin outcome
###############
# BUG A. This study was run with the CONTINUOUS coefficient table
# (b0 = c(0.4, 0.2, ...), b1 = -0.05, b2 = c(2, 2, ...), plus the s2/s4/s5
# variance parameters) fed through plogis to produce a binary outcome. Its
# results are therefore not comparable with the main binary study, and with
# b2 = 2 on a logit scale the outcome model saturates.
#
# set = "binary_ci" resolves to the continuous table while LEGACY_BIN_CI_PARAMS
# is TRUE in R/dgm_scenarios.R, and to the binary table once it is FALSE.
# Flipping it invalidates everything under ../results/confidence_intervals/binary
# and requires the study (and optimal_sf/bin_ci_sf) to re-run.

source(here::here("R", "dgm_scenarios.R"))

generate_binary_scenario_data <- function(scenario, n, return_truth = TRUE,
                                          seed = NULL) {
  generate_scenario_data(scenario, n, set = "binary_ci",
                         return_truth = return_truth, seed = seed)
}

get_binary_oracle_info <- function(scenario, bW) {
  get_oracle_info(scenario, bW, set = "binary_ci")
}
