##########
# title: correlated binary study - the one definition of its parameter grid
##########
# Sourced by bin_corr_analysis.R, bin_corr_check.R, bin_corr_collect.R and
# bin_corr_metrics.R. The array index handed to bin_corr_analysis.R is a row
# number of `grid`, so `grid` must never be filtered or reordered after
# construction.
#
# binary/'s scenarios 1-4 with correlated covariates, at every rho in
# CORR_RHOS (R/dgm_scenarios.R). rho varies slowest, so rows 1-8000 are rho = 0
# (jobscripts/bin_corr_1.sh) and rows 8001-16000 rho = 0.5 (bin_corr_2.sh). A
# run is seeded by its run index alone, so the same run at the two rhos is a
# paired dataset - see ../README.md.

library(here)
source(here("R", "pipeline.R"))
source(here("R", "dgm_scenarios.R"))

study <- study_config(
  name     = "correlated/binary",
  prefix   = "bin_corr",
  res_path = file.path(dirname(here()), "results", "correlated", "binary"),
  grid = expand.grid(
    scenario = c(1:4),
    n = c(100, 250, 500, 1000),
    run = c(1:500),
    rho = CORR_RHOS,
    stringsAsFactors = FALSE
  ),
  path_cols   = c("rho", "scenario", "n"),
  path_prefix = c(rho = "rho_", scenario = "scenario_"),
  n_sims      = 500,
  failed_file = here("sample_size", "correlated", "binary", "jobscripts",
                     "failed_ids.txt")
)
