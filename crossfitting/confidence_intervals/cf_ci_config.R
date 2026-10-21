##########
# title: crossfitting CI pilot - the one definition of its parameter grid
##########
# Sourced by study.R check. Kept separate from crossfitting/cf_config.R, with
# its own jobscripts/failed_ids.txt, so a re-run of one study never clobbers the
# other's todo list.
#
# The grid here must stay byte-for-byte identical to the one inline in
# cf_ci_analysis.R, which qsub calls on the cluster and keeps its own copy.

library(here)
source(here("R", "pipeline.R"))

study <- study_config(
  name     = "crossfitting/confidence_intervals",
  prefix   = "cf_ci",
  res_path = file.path(dirname(here()), "results", "crossfitting_ci"),
  grid = expand.grid(
    scenario = c(1, 4, 8),
    n = 500,
    run = c(1:50),
    stringsAsFactors = FALSE
  ),
  path_cols   = "scenario",
  n_sims      = 50,
  failed_file = here("crossfitting", "confidence_intervals", "jobscripts",
                      "failed_ids.txt")
)
