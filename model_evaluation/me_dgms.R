###############
# title: data generating process - model evaluation study
###############
# The scenario tables and the generator live in R/dgm_scenarios.R. This file
# names this study's slice of them - the continuous_corr_0.5 set: continuous/'s 10-scenario
# table with the covariates drawn from the copula at latent correlation
# ME_RHO, the same set sample_size/ uses as its primary arm. tau(x)
# and m0(x) are unchanged; only the covariates' joint distribution moves, and
# with it bW, the ATE and SD(tau) wherever g multiplies two modifiers
# (scenario 3 here, X4*X5) - see sample_size/README.md. This study
# only ever runs 4 of the scenarios, see me_config.R.
#
# Only rho = 0.5 is run, so it is a constant here rather than a grid column:
# the grid, the array indices and the result paths keep their shape.
#
# No get_*_oracle_info() wrapper, unlike sample_size/'s _dgms.R files: none of
# the 9 candidates is an oracle.

source(here::here("R", "dgm_scenarios.R"))

ME_RHO <- 0.5

# calib_n: the n bW is calibrated at (see generate_scenario_data()). The study
# leaves it at n; me_infl_investigation.R draws a population-size sample but
# keeps a study n here, so the true ATE is the study's. Consumes no RNG.
generate_me_scenario_data <- function(scenario, n, return_truth = TRUE, seed = NULL,
                                      calib_n = n) {
  generate_scenario_data(scenario, n, set = corr_set("continuous", ME_RHO),
                         return_truth = return_truth, seed = seed, calib_n = calib_n)
}
