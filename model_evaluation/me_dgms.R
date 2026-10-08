###############
# title: data generating process - model evaluation study
###############
# The scenario tables and the generator live in R/dgm_scenarios.R. This file
# names this study's slice of them - the continuous_corr_0.5 set (since
# 2026-10-08; it ran on set = "continuous" before): continuous/'s 10-scenario
# table with the covariates drawn from the copula at latent correlation
# ME_RHO, the same set sample_size/correlated/ uses as its primary arm. tau(x)
# and m0(x) are unchanged; only the covariates' joint distribution moves, and
# with it bW, the ATE and SD(tau) wherever g multiplies two modifiers
# (scenario 8 here, X3*X4) - see sample_size/correlated/README.md. This study
# only ever runs 4 of the scenarios, see me_config.R.
#
# Only rho = 0.5 is run, so it is a constant here rather than a grid column:
# the grid, the array indices and the result paths keep their shape.
#
# Replaces the benchtm::generate_scen_data() call this study used to make:
# gen$truth$tau (a data.frame(p0, p1, tau)) is the direct replacement for
# benchtm's adat$trt_effect, and generate_scenario_data() already returns Y,
# W as the first two columns, so results$truth <- gen$truth saves verbatim.
#
# No get_*_oracle_info() wrapper here, unlike sample_size/{continuous,binary}'s _dgms.R
# files - none of this study's 9 candidate models are oracle estimators, so
# a wrapper nothing calls would be dead code from the moment it's written.

source(here::here("R", "dgm_scenarios.R"))

ME_RHO <- 0.5

generate_me_scenario_data <- function(scenario, n, return_truth = TRUE, seed = NULL) {
  generate_scenario_data(scenario, n, set = corr_set("continuous", ME_RHO),
                         return_truth = return_truth, seed = seed)
}
