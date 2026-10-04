##########
# title: correlated binary study - the one definition of its parameter grid
##########
# Sourced by bin_corr_analysis.R, bin_corr_check.R, bin_corr_collect.R and
# bin_corr_metrics.R. The array index handed to bin_corr_analysis.R is a row
# number of `grid`, so `grid` must never be filtered or reordered after
# construction.
#
# binary/'s scenarios 1-10 with correlated covariates, at every rho in
# CORR_RHOS (R/dgm_scenarios.R). A run is seeded by its run index alone, so the
# same run at the two rhos is a paired dataset - see ../README.md.
#
# Scenarios 1-4 run to 500 runs, rows 1-16000. rho varies slowest, so rows
# 1-8000 are rho = 0 (jobscripts/bin_corr_1.sh) and rows 8001-16000 rho = 0.5
# (bin_corr_2.sh). Scenarios 5-10 (the appendix's, since 2026-10-04) run to
# 100, like binary/'s, and are appended as a second block rather than widening
# `scenario` above, so rows 1-16000 keep the meaning they were submitted under
# and the new rows are one contiguous range, 16001-20800, for
# jobscripts/bin_corr_extra.sh: rows 16001-18400 are rho = 0, 18401-20800
# rho = 0.5.

library(here)
source(here("R", "pipeline.R"))
source(here("R", "dgm_scenarios.R"))

study <- study_config(
  name     = "correlated/binary",
  prefix   = "bin_corr",
  res_path = file.path(dirname(here()), "results", "correlated", "binary"),
  grid = rbind(
    expand.grid(
      scenario = c(1:4),
      n = c(100, 250, 500, 1000),
      run = c(1:500),
      rho = CORR_RHOS,
      stringsAsFactors = FALSE
    ),
    expand.grid(
      scenario = c(5:10),
      n = c(100, 250, 500, 1000),
      run = c(1:100),
      rho = CORR_RHOS,
      stringsAsFactors = FALSE
    )
  ),
  path_cols   = c("rho", "scenario", "n"),
  path_prefix = c(rho = "rho_", scenario = "scenario_"),
  n_sims      = 500,
  failed_file = here("sample_size", "correlated", "binary", "jobscripts",
                     "failed_ids.txt")
)
