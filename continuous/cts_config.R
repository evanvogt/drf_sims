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
  prefix   = "cts",
  res_path = file.path(dirname(here()), "results", "continuous"),
  # Scenarios 1-4 go to 500 runs. Their runs 101-500 are appended as
  # a second block rather than widening `run` above, so rows 1-4000 keep the
  # meaning they were submitted under and the new rows are one contiguous range,
  # 4001-10400, for jobscripts/cts_extra.sh.
  grid = rbind(
    expand.grid(
      scenario = c(1:10),
      n = c(100, 250, 500, 1000),
      run = c(1:100),
      stringsAsFactors = FALSE
    ),
    expand.grid(
      scenario = c(1:4),
      n = c(100, 250, 500, 1000),
      run = c(101:500),
      stringsAsFactors = FALSE
    )
  ),
  path_cols   = c("scenario", "n"),
  n_sims      = 100,
  failed_file = here("continuous", "jobscripts", "failed_ids.txt")
)
