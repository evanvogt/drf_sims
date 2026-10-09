##########
# title: interim-analysis validation study - the one definition of its parameter grid
##########
# Sourced by cts_val_analysis.R, cts_val_check.R, cts_val_collect.R and
# cts_val_metrics.R. The array index handed to cts_val_analysis.R is a row
# number of `grid`, so `grid` must never be filtered or reordered after
# construction.
#
# Correlated covariates only, at rho = 0.5 - the correlated studies' primary arm
# (sample_size/correlated/, R/dgm_scenarios.R's CORRELATED COVARIATES). rho is a
# path column with the same "rho_" prefix those studies use, so results sit
# under rho_0.5/ and a rho = 0 arm would be a one-value change to the grid.

library(here)
source(here("R", "pipeline.R"))

study <- study_config(
  name     = "validation/continuous",
  prefix   = "cts_val",
  res_path = file.path(dirname(here()), "results", "validation", "continuous"),
  grid = expand.grid(
    rho = c(0, 0.5),
    scenario = 2,
    n = 1000,
    # round() is load-bearing, not cosmetic: interim_prop is a path_cols entry,
    # so its as.character() form becomes a results directory name. Unrounded,
    # seq(by = 0.05) yields values like 0.30000000000000004 and get_results()
    # can no longer match the directory back to its grid row.
    interim_prop = round(seq(0.25, 0.75, by = 0.05), 2),
    run = c(1:100),
    stringsAsFactors = FALSE
  ),
  path_cols   = c("rho", "scenario", "n", "interim_prop"),
  path_prefix = c(rho = "rho_", scenario = "scenario_"),
  n_sims      = 100,
  failed_file = here("validation", "continuous", "jobscripts", "failed_ids.txt")
)
