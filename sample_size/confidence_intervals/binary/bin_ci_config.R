##########
# title: binary CI study - the one definition of its parameter grid
##########
# NOTE: bug A (continuous coefficient table) was fixed; study now always uses binary table.

library(here)
source(here("R", "pipeline.R"))

study <- study_config(
  name     = "confidence_intervals/binary",
  prefix   = "bin_ci",
  res_path = file.path(dirname(here()), "results", "confidence_intervals", "binary"),
  # 100 runs for every scenario. Scenarios 1-4 used to go to 500 (rows
  # 20001-52000, jobscripts/bin_ci_extra_{1..4}.sh); cut back to 100 on
  # 2026-09-27, so rows 1-20000 mean what they always did.
  grid = expand.grid(
    scenario = c(1:10),
    n = c(500, 1000),
    CI_sf = seq(0.05, 0.5, 0.05),
    run = c(1:100),
    stringsAsFactors = FALSE
  ),
  path_cols   = c("scenario", "n", "CI_sf"),
  n_sims      = 100,
  failed_file = here("sample_size", "confidence_intervals", "binary", "jobscripts", "failed_ids.txt")
)
