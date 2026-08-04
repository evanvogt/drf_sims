##########
# title: continuous study - the one definition of its parameter grid
##########
# Sourced by cts_analysis.R, cts_check.R, cts_collect.R and cts_metrics.R.
# The array index handed to cts_analysis.R is a row number of `grid`, so `grid`
# must never be filtered or reordered after construction.

library(here)
source(here("R", "pipeline.R"))

study <- study_config(
  name     = "continuous",
  res_path = file.path(dirname(here()), "results", "continuous"),
  grid = expand.grid(
    scenario = c(1:10),
    n = c(100, 250, 500, 1000),
    run = c(1:100),
    stringsAsFactors = FALSE
  ),
  path_cols   = c("scenario", "n"),
  n_sims      = 100,
  failed_file = here("continuous", "jobscripts", "failed_ids.txt")
)
