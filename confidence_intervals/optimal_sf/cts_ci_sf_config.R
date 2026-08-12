##########
# title: optimal sample.fraction calibration (continuous) - the one definition of its parameter grid
##########
# Sourced by cts_ci_sf_check.R. The grid here must stay byte-for-byte
# identical to the one inline in cts_ci_sf_analysis.R - that script is what
# qsub calls on the cluster and is deliberately left untouched.

library(here)
source(here("R", "pipeline.R"))

study <- study_config(
  name     = "confidence_intervals/optimal_sf (cts)",
  res_path = file.path(dirname(here()), "results", "confidence_intervals",
                        "continuous", "sf_calibration"),
  grid = expand.grid(
    scenario = c(1:10),
    n        = c(500, 1000),
    run      = c(1:100),
    stringsAsFactors = FALSE
  ),
  path_cols   = c("scenario", "n"),
  n_sims      = 100,
  failed_file = here("confidence_intervals", "optimal_sf", "jobscripts",
                      "failed_cts_ids.txt")
)
