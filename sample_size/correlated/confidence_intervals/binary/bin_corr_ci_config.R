##########
# title: correlated binary CI study - the one definition of its parameter grid
##########
# Sourced by bin_corr_ci_analysis.R, bin_corr_ci_check.R, bin_corr_ci_collect.R
# and bin_corr_ci_metrics.R. The array index handed to bin_corr_ci_analysis.R is
# a row number of `grid`, so `grid` must never be filtered or reordered after
# construction.
#
# confidence_intervals/binary/'s design (the full CI_sf sweep, n in
# {500, 1000}, 100 runs) on correlated/'s scenarios 1-4, at every rho in
# CORR_RHOS (R/dgm_scenarios.R). rho varies slowest, so rows 1-8000 are rho = 0
# and rows 8001-16000 rho = 0.5. The PBS array cap is 10,000, so the jobscripts
# split at 10000/10001 instead (bin_corr_ci_1.sh, bin_corr_ci_2.sh). A run is
# seeded by its run index alone, so the same run at the two rhos - and at every
# CI_sf - is a paired dataset; see ../README.md.

library(here)
source(here("R", "pipeline.R"))
source(here("R", "dgm_scenarios.R"))

study <- study_config(
  name     = "correlated/confidence_intervals/binary",
  prefix   = "bin_corr_ci",
  res_path = file.path(dirname(here()), "results", "correlated",
                       "confidence_intervals", "binary"),
  grid = expand.grid(
    scenario = c(1:4),
    n = c(500, 1000),
    CI_sf = seq(0.05, 0.5, 0.05),
    run = c(1:100),
    rho = CORR_RHOS,
    stringsAsFactors = FALSE
  ),
  path_cols   = c("rho", "scenario", "n", "CI_sf"),
  path_prefix = c(rho = "rho_", scenario = "scenario_"),
  n_sims      = 100,
  failed_file = here("sample_size", "correlated", "confidence_intervals",
                     "binary", "jobscripts", "failed_ids.txt")
)
