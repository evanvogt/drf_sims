##########
# title: competing risks study - the one definition of its parameter grid
##########
# Results are laid out as scenario_<k>/<n>/censor_<TRUE|FALSE>, hence the
# path_prefix entry for `censoring`.

library(here)
source(here("R", "pipeline.R"))

study <- study_config(
  name     = "competing_risk",
  prefix   = "surv",
  res_path = file.path(dirname(here()), "results", "competing_risk"),
  grid = expand.grid(
    scenario = 1:7,
    censoring = c(TRUE, FALSE),
    n = c(500),
    run = 1:100,
    stringsAsFactors = FALSE
  ),
  path_cols   = c("scenario", "n", "censoring"),
  path_prefix = c(scenario = "scenario_", censoring = "censor_"),
  n_sims      = 100,
  failed_file = here("competing_risk", "jobscripts", "failed_ids.txt")
)
