##########
# title: missing-covariate binary study - the one parameter grid
##########
# As missing/continuous, and with the same bugs C and D fixed by construction
# (bug B was continuous-only - this study already spelled the mechanisms MNAR).
#
# NOTE: this study's DGM carried three further defects (continuous coefficients,
# continuous power calibration, log-odds truth), all fixed together.
# See missing/binary/README.md.

library(here)
source(here("R", "pipeline.R"))

MISS_METHODS_STUDY <- c("complete_cases", "mean_imputation", "missforest",
                        "regression", "missing_indicator", "IPW",
                        "multiple_imputation", "none", "complete_data")

grid <- expand.grid(
  scenario  = c(1, 3, 4, 5),
  n         = c(500),
  type      = c("both"),
  prop      = c(0.3),
  mechanism = c("MAR", "MNAR", "MNAR-Y"),
  method    = MISS_METHODS_STUDY,
  run       = c(1:100),
  stringsAsFactors = FALSE
)

grid <- grid[!(grid$scenario == 1 & grid$mechanism == "MNAR-Y"), ]

# Scenario 2 appended as a second block, rows 9901-12600, for
# jobscripts/bin_miss_extra.sh - see missing/continuous/cts_miss_config.R,
# including why rows 1-9900 changed meaning at the 2026-09-26 renumbering.
grid <- rbind(grid, expand.grid(
  scenario  = 2,
  n         = c(500),
  type      = c("both"),
  prop      = c(0.3),
  mechanism = c("MAR", "MNAR", "MNAR-Y"),
  method    = MISS_METHODS_STUDY,
  run       = c(1:100),
  stringsAsFactors = FALSE
))
rownames(grid) <- NULL

study <- study_config(
  name     = "missing/binary",
  prefix   = "bin_miss",
  res_path = file.path(dirname(here()), "results", "missing", "binary"),
  grid        = grid,
  path_cols   = c("scenario", "n", "type", "prop", "mechanism", "method"),
  n_sims      = 100,
  failed_file = here("missing", "binary", "jobscripts", "failed_ids.txt")
)
