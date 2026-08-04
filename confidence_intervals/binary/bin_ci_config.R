##########
# title: binary CI study - the one definition of its parameter grid
##########
# NOTE: this study's DGM carries bug A - see confidence_intervals/binary/README.md
# and LEGACY_BIN_CI_PARAMS in R/dgm_scenarios.R. Fixing it invalidates everything
# under ../results/confidence_intervals/binary.

library(here)
source(here("R", "pipeline.R"))

study <- study_config(
  name     = "confidence_intervals/binary",
  res_path = file.path(dirname(here()), "results", "confidence_intervals", "binary"),
  grid = expand.grid(
    scenario = c(1:10),
    n = c(500, 1000),
    CI_sf = seq(0.05, 0.5, 0.05),
    run = c(1:100),
    stringsAsFactors = FALSE
  ),
  path_cols   = c("scenario", "n", "CI_sf"),
  n_sims      = 100,
  failed_file = here("confidence_intervals", "binary", "jobscripts", "failed_ids.txt")
)
