###############
# title: data generating process - correlated binary study
###############
# The scenario tables and the generator live in R/dgm_scenarios.R. This file
# names this study's slice of them: the binary_corr_<rho> sets, binary/'s
# scenarios 1-4 (risk-difference scale, RD_SCALE unchanged) with the covariates
# drawn from the copula at latent correlation rho. The oracle formula returns
# the risk itself.

source(here::here("R", "dgm_scenarios.R"))

generate_binary_scenario_data <- function(scenario, n, rho,
                                          return_truth = TRUE, seed = NULL) {
  generate_scenario_data(scenario, n, set = corr_set("binary", rho),
                         return_truth = return_truth, seed = seed)
}

get_binary_oracle_info <- function(scenario, bW, rho) {
  get_oracle_info(scenario, bW, set = corr_set("binary", rho))
}
