##########
# title: missing-covariate continuous study - the one parameter grid
##########
# Every script reads the grid from here, which rules out bugs B-D
# (docs/DEVLOG.md):
#   - one spelling of each mechanism (R/dgm_scenarios.R's MISS_MECHS)
#   - "complete_data" stays in the grid: it is the reference arm, with no
#     missingness introduced, that the other methods are scored against
#   - the grid is never filtered after construction; to run a subset, select
#     indices with grid_indices() so row numbers keep their meaning

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
  mechanism = c("MAR", "MNAR-Y0", "MNAR-tau"),
  method    = MISS_METHODS_STUDY,
  run       = c(1:100),
  stringsAsFactors = FALSE
)

# scenario 1 has no HTE, so a mechanism that depends on the treatment effect
# through U is not defined for it. Applied at construction, before any row
# numbers are handed out.
grid <- grid[!(grid$scenario == 1 & grid$mechanism == "MNAR-tau"), ]

# Scenario 2 is appended as its own block, not added to `scenario` above, so its
# rows are one contiguous range, 9901-12600, for jobscripts/cts_miss_extra.sh.
grid <- rbind(grid, expand.grid(
  scenario  = 2,
  n         = c(500),
  type      = c("both"),
  prop      = c(0.3),
  mechanism = c("MAR", "MNAR-Y0", "MNAR-tau"),
  method    = MISS_METHODS_STUDY,
  run       = c(1:100),
  stringsAsFactors = FALSE
))
rownames(grid) <- NULL

study <- study_config(
  name     = "missing/continuous",
  prefix   = "cts_miss",
  res_path = file.path(dirname(here()), "results", "missing", "continuous"),
  grid        = grid,
  path_cols   = c("scenario", "n", "type", "prop", "mechanism", "method"),
  n_sims      = 100,
  failed_file = here("missing", "continuous", "jobscripts", "failed_ids.txt")
)
