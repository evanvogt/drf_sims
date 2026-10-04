##########
# title: single-event survival study - the one definition of its parameter grid
##########
# Sourced by se_analysis.R, se_check.R, se_collect.R and se_metrics.R. The
# array index handed to se_analysis.R is a row number of `grid`, so `grid` must
# never be filtered or reordered after construction.
#
# Scenario varies fastest, then censoring: index 3 is scenario 3 with censoring
# on, run 1, and index 6 the same run with censoring off.
#
# Results live in results/single_event/, NOT under results/competing_risk/: the
# parent README's archive step tars and removes that whole tree.

library(here)
source(here("R", "pipeline.R"))

study <- study_config(
  name     = "competing_risk/single_event",
  prefix   = "se",
  res_path = file.path(dirname(here()), "results", "single_event"),
  grid = expand.grid(
    scenario = 1:3,
    censoring = c(TRUE, FALSE),
    n = c(500),
    run = 1:500,
    stringsAsFactors = FALSE
  ),
  path_cols   = c("scenario", "n", "censoring"),
  path_prefix = c(scenario = "scenario_", censoring = "censor_"),
  n_sims      = 500,
  failed_file = here("competing_risk", "single_event", "jobscripts",
                     "failed_ids.txt")
)
