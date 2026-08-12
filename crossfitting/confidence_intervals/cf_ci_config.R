##########
# title: crossfitting CI pilot - the one definition of its parameter grid
##########
# Sourced by cf_ci_check.R. Kept separate from crossfitting/cf_config.R (a
# different failed_file - failed_ci_ids.txt) so a re-run of one study never
# clobbers the other's todo list, matching the previous check script's intent.
#
# The grid here must stay byte-for-byte identical to the one inline in
# cf_ci_analysis.R - that script is what qsub calls on the cluster and is
# deliberately left untouched.

library(here)
source(here("R", "pipeline.R"))

study <- study_config(
  name     = "crossfitting/confidence_intervals",
  res_path = file.path(dirname(here()), "results", "crossfitting_ci"),
  grid = expand.grid(
    scenario = c(1, 6, 9),
    n = 500,
    run = c(1:50),
    stringsAsFactors = FALSE
  ),
  path_cols   = c("scenario", "n"),
  n_sims      = 50,
  failed_file = here("crossfitting", "confidence_intervals", "jobscripts",
                      "failed_ci_ids.txt")
)
