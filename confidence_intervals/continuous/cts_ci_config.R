##########
# title: continuous CI study - the one definition of its parameter grid
##########
# CI_sf is swept because the half-sample bootstrap's sample.fraction has to be
# calibrated; confidence_intervals/optimal_sf/ picks the value from this sweep.

library(here)
source(here("R", "pipeline.R"))

study <- study_config(
  name     = "confidence_intervals/continuous",
  res_path = file.path(dirname(here()), "results", "confidence_intervals", "continuous"),
  grid = expand.grid(
    scenario = c(1:10),
    n = c(500, 1000),
    CI_sf = seq(0.05, 0.5, 0.05),
    run = c(1:100),
    stringsAsFactors = FALSE
  ),
  path_cols   = c("scenario", "n", "CI_sf"),
  n_sims      = 100,
  failed_file = here("confidence_intervals", "continuous", "jobscripts", "failed_ids.txt")
)
