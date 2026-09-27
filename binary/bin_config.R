##########
# title: binary study - the one definition of its parameter grid
##########
# Sourced by bin_analysis.R, bin_check.R, bin_collect.R and bin_metrics.R.
#
# GRID DISAGREEMENT, since resolved. Before this file existed the grid was
# declared three times and they did not agree:
#
#   bin_analysis.R          scenario = c(1:10)      -> 4000 rows
#   bin_check/bin_collect   scenario = c(1, 3, 8, 9) -> 1600 rows
#   jobscripts/bin_1.sh     #PBS -J 1-1600
#
# 4 scenarios x 4 sample sizes x 100 runs is exactly 1600, so c(1, 3, 8, 9) was
# the intent and the analysis script's c(1:10) was stale. Because expand.grid
# varies the first column fastest, submitting indices 1-1600 against the
# 4000-row grid actually ran runs 1-40 of ALL TEN scenarios. bin_collect.R then
# looked for scenarios 1, 3, 8 and 9 and found 40 runs in each, so the study
# quietly has 40 replicates per cell rather than the intended 100.
#
# Anything already under ../results/binary was produced by the old mapping.
# The study re-runs anyway (bug F, bug P), on the grid below: all ten
# scenarios at 100 runs, with 1-4 taken to 500. (Scenario numbers above are
# the pre-2026-09-26 ones; old 1, 3, 8, 9 are now 1-4 - see R/dgm_scenarios.R.)

library(here)
source(here("R", "pipeline.R"))

study <- study_config(
  name     = "binary",
  prefix   = "bin",
  res_path = file.path(dirname(here()), "results", "binary"),
  # Scenarios 1-4 go to 500 runs. Their runs 101-500 are appended as
  # a second block rather than widening `run` above, so rows 1-4000 keep the
  # meaning they were submitted under and the new rows are one contiguous range,
  # 4001-10400, for jobscripts/bin_extra.sh.
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
  failed_file = here("binary", "jobscripts", "failed_ids.txt")
)
