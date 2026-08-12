##########
# title: crossfitting comparison study - the one definition of its parameter grid
##########
# Sourced by cf_check.R (and, indirectly, R/study_registry.R). The grid here
# must stay byte-for-byte identical to the one inline in cf_analysis.R (same
# columns, same order, same values) - cf_analysis.R is what qsub calls on the
# cluster and is deliberately left untouched, so the array index it uses and
# the row number this file hands out must keep meaning the same thing.

library(here)
source(here("R", "pipeline.R"))

study <- study_config(
  name     = "crossfitting",
  prefix   = "cf",
  res_path = file.path(dirname(here()), "results", "crossfitting"),
  grid = expand.grid(
    scenario = c(1, 4, 6, 9),
    n = 500,
    run = c(1:500),
    stringsAsFactors = FALSE
  ),
  path_cols   = c("scenario", "n"),
  n_sims      = 500,
  failed_file = here("crossfitting", "jobscripts", "failed_ids.txt")
)
