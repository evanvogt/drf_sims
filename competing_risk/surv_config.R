##########
# title: competing risks study - the one definition of its parameter grid
##########
# Results are laid out as rho_<rho>/scenario_<k>/<n>/censor_<TRUE|FALSE>, hence
# the path_prefix entries for `rho` and `censoring`.
#
# Covariates are drawn from the copula at every rho in CORR_RHOS
# (R/dgm_scenarios.R; see surv_dgm.R). rho varies slowest, so rows 1-7000 are
# rho = 0, the primary analysis (jobscripts/surv_1.sh), and rows 7001-14000 are
# rho = 0.5 (surv_2.sh). A run is seeded by its run index alone, so the same
# run at the two rhos is a paired dataset. The array index is a row number of
# `grid`, so `grid` must never be filtered or reordered after construction.

library(here)
source(here("R", "pipeline.R"))
source(here("R", "dgm_scenarios.R"))

study <- study_config(
  name     = "competing_risk",
  prefix   = "surv",
  res_path = file.path(dirname(here()), "results", "competing_risk"),
  grid = expand.grid(
    scenario = 1:7,
    censoring = c(TRUE, FALSE),
    n = c(500),
    run = 1:500,
    rho = CORR_RHOS,
    stringsAsFactors = FALSE
  ),
  path_cols   = c("rho", "scenario", "n", "censoring"),
  path_prefix = c(rho = "rho_", scenario = "scenario_", censoring = "censor_"),
  n_sims      = 500,
  failed_file = here("competing_risk", "jobscripts", "failed_ids.txt")
)
