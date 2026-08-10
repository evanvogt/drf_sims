###############
# title: data generating process - model evaluation study
###############
# The scenario tables and the generator live in R/dgm_scenarios.R. This file
# names this study's slice of them - set = "continuous", the same 10-scenario
# table continuous/ uses (this study only ever runs 4 of them, see
# me_config.R, but the underlying scenario table is unchanged).
#
# Replaces the benchtm::generate_scen_data() call this study used to make:
# gen$truth$tau (a data.frame(p0, p1, tau)) is the direct replacement for
# benchtm's adat$trt_effect, and generate_scenario_data() already returns Y,
# W as the first two columns, so results$truth <- gen$truth saves verbatim.
#
# No get_*_oracle_info() wrapper here, unlike continuous/binary's _dgms.R
# files - none of this study's 9 candidate models are oracle estimators, so
# a wrapper nothing calls would be dead code from the moment it's written.

source(here::here("R", "dgm_scenarios.R"))

generate_me_scenario_data <- function(scenario, n, return_truth = TRUE, seed = NULL) {
  generate_scenario_data(scenario, n, set = "continuous",
                         return_truth = return_truth, seed = seed)
}
