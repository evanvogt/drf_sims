###############
# title: data generating process - correlated continuous study
###############
# The scenario tables and the generator live in R/dgm_scenarios.R. This file
# names this study's slice of them: the continuous_corr_<rho> sets, continuous/'s
# scenarios 1-4 with the covariates drawn from the copula at latent
# correlation rho.

source(here::here("R", "dgm_scenarios.R"))

generate_continuous_scenario_data <- function(scenario, n, rho,
                                              return_truth = TRUE, seed = NULL) {
  generate_scenario_data(scenario, n, set = corr_set("continuous", rho),
                         return_truth = return_truth, seed = seed)
}

get_continuous_oracle_info <- function(scenario, bW, rho) {
  get_oracle_info(scenario, bW, set = corr_set("continuous", rho))
}
