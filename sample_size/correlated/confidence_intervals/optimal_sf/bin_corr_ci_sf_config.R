##########
# title: correlated optimal sample.fraction calibration (binary) - the one definition of its parameter grid
##########
# Sourced by bin_corr_ci_sf_analysis.R, bin_corr_ci_sf_check.R,
# bin_corr_ci_sf_collect.R and bin_corr_ci_sf_metrics.R. The array index is a
# row number of `grid`, so `grid` must never be filtered or reordered.
#
# confidence_intervals/optimal_sf/'s design (scenarios 1-4, n in {500, 1000},
# 100 runs, no CI_sf axis - that is what each run picks) on correlated/'s
# copula, at every rho in CORR_RHOS. rho varies slowest: rows 1-800 are
# rho = 0, rows 801-1600 rho = 0.5, and the same run at the two rhos is a
# paired dataset (../README.md).
#
# failed_bin_ids.txt, not failed_ids.txt: jobscripts/ serves both optimal_sf
# studies, as in the parent, so the two todo lists are named apart.

library(here)
source(here("R", "pipeline.R"))
source(here("R", "dgm_scenarios.R"))

study <- study_config(
  name     = "correlated/confidence_intervals/optimal_sf (bin)",
  prefix   = "bin_corr_ci_sf",
  res_path = file.path(dirname(here()), "results", "correlated",
                       "confidence_intervals", "binary", "sf_calibration"),
  grid = expand.grid(
    scenario = c(1:4),
    n        = c(500, 1000),
    run      = c(1:100),
    rho      = CORR_RHOS,
    stringsAsFactors = FALSE
  ),
  path_cols   = c("rho", "scenario", "n"),
  path_prefix = c(rho = "rho_", scenario = "scenario_"),
  n_sims      = 100,
  failed_file = here("sample_size", "correlated", "confidence_intervals",
                     "optimal_sf", "jobscripts", "failed_bin_ids.txt")
)
