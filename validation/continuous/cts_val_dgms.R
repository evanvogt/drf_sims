###############
# title: data generating process - interim-analysis validation study
###############
# The scenario tables and the generator live in R/dgm_scenarios.R. This file
# names this study's slice of them: the continuous_corr_<rho> sets that
# sample_size/correlated/continuous/cts_corr_dgms.R also names - continuous/'s
# scenarios with X1-X5 and X01-X03 drawn from a Gaussian copula with
# exchangeable latent correlation rho. It replaces the standalone fork that
# used to live in cts_dgm_validation.R.

source(here::here("R", "dgm_scenarios.R"))

generate_continuous_scenario_data <- function(scenario, n, rho,
                                              return_truth = TRUE, seed = NULL) {
  generate_scenario_data(scenario, n, set = corr_set("continuous", rho),
                         return_truth = return_truth, seed = seed)
}

get_continuous_oracle_info <- function(scenario, bW, rho) {
  get_oracle_info(scenario, bW, set = corr_set("continuous", rho))
}
