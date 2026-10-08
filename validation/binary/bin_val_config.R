##########
# title: interim-analysis validation study, binary outcome - the one definition
#        of its parameter grid
##########
# Sourced by bin_val_analysis.R, bin_val_check.R, bin_val_collect.R and
# bin_val_metrics.R. The array index handed to bin_val_analysis.R is a row
# number of `grid`, so `grid` must never be filtered or reordered after
# construction.
#
# The same grid as continuous/cts_val_config.R, on the binary_corr_0.5 set:
# correlated covariates at rho = 0.5, scenario 2, one trial of 1000 split at
# eleven interim points, 100 runs. rho is a path column with the "rho_" prefix
# the correlated studies use, so results sit under rho_0.5/.

library(here)
source(here("R", "pipeline.R"))

study <- study_config(
  name     = "validation/binary",
  prefix   = "bin_val",
  res_path = file.path(dirname(here()), "results", "validation", "binary"),
  grid = expand.grid(
    rho = 0.5,
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
  failed_file = here("validation", "binary", "jobscripts", "failed_ids.txt")
)
