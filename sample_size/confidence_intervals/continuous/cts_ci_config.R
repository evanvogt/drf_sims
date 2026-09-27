##########
# title: continuous CI study - the one definition of its parameter grid
##########
# CI_sf is swept because the half-sample bootstrap's sample.fraction has to be
# calibrated; confidence_intervals/optimal_sf/ picks the value from this sweep.

library(here)
source(here("R", "pipeline.R"))

study <- study_config(
  name     = "confidence_intervals/continuous",
  prefix   = "cts_ci",
  res_path = file.path(dirname(here()), "results", "confidence_intervals", "continuous"),
  # Scenarios 1-4 go to 500 runs. Their runs 101-500 are appended as
  # a second block rather than widening `run` above, so rows 1-20000 keep the
  # meaning they were submitted under and the new rows are one contiguous range,
  # 20001-52000, for jobscripts/cts_ci_extra_{1..4}.sh.
  grid = rbind(
    expand.grid(
      scenario = c(1:10),
      n = c(500, 1000),
      CI_sf = seq(0.05, 0.5, 0.05),
      run = c(1:100),
      stringsAsFactors = FALSE
    ),
    expand.grid(
      scenario = c(1:4),
      n = c(500, 1000),
      CI_sf = seq(0.05, 0.5, 0.05),
      run = c(101:500),
      stringsAsFactors = FALSE
    )
  ),
  path_cols   = c("scenario", "n", "CI_sf"),
  n_sims      = 100,
  failed_file = here("sample_size", "confidence_intervals", "continuous", "jobscripts", "failed_ids.txt")
)
