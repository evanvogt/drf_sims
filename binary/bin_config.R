##########
# title: binary study - the one definition of its parameter grid
##########
# Sourced by bin_analysis.R, bin_check.R, bin_collect.R and bin_metrics.R.
#
# GRID DISAGREEMENT, now resolved in favour of c(1, 3, 8, 9). Before this file
# existed the grid was declared three times and they did not agree:
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
# The study re-runs anyway for bug F, and this grid gives the intended design.

library(here)
source(here("R", "pipeline.R"))

study <- study_config(
  name     = "binary",
  res_path = file.path(dirname(here()), "results", "binary"),
  grid = expand.grid(
    scenario = c(1, 3, 8, 9),
    n = c(100, 250, 500, 1000),
    run = c(1:100),
    stringsAsFactors = FALSE
  ),
  path_cols   = c("scenario", "n"),
  n_sims      = 100,
  failed_file = here("binary", "jobscripts", "failed_ids.txt")
)
