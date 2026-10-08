###############
# title: data generating process - interim-analysis validation study, binary
###############
# The scenario tables and the generator live in R/dgm_scenarios.R. This file
# names this study's slice of them: the binary_corr_<rho> sets that
# sample_size/correlated/binary/bin_corr_dgms.R also names - binary/'s
# scenarios on the risk-difference scale, with X1-X5 and X01-X03 drawn from a
# Gaussian copula with exchangeable latent correlation rho.

source(here::here("R", "dgm_scenarios.R"))

generate_binary_scenario_data <- function(scenario, n, rho,
                                          return_truth = TRUE, seed = NULL) {
  generate_scenario_data(scenario, n, set = corr_set("binary", rho),
                         return_truth = return_truth, seed = seed)
}

get_binary_oracle_info <- function(scenario, bW, rho) {
  get_oracle_info(scenario, bW, set = corr_set("binary", rho))
}
