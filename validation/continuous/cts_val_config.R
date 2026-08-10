##########
# title: interim-analysis validation study - the one definition of its parameter grid
##########
# Sourced by cts_val_analysis.R, cts_val_check.R, cts_val_collect.R and
# cts_val_metrics.R. The array index handed to cts_val_analysis.R is a row
# number of `grid`, so `grid` must never be filtered or reordered after
# construction.

library(here)
source(here("R", "pipeline.R"))

study <- study_config(
  name     = "validation/continuous",
  res_path = file.path(dirname(here()), "results", "validation", "continuous"),
  grid = expand.grid(
    scenario = 3,
    n = 1000,
    interim_prop = c(0.25, 0.5, 0.75),
    run = c(1:100),
    stringsAsFactors = FALSE
  ),
  path_cols   = c("scenario", "n", "interim_prop"),
  n_sims      = 100,
  failed_file = here("validation", "continuous", "jobscripts", "failed_ids.txt")
)
