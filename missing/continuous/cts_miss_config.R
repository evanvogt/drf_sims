##########
# title: missing-covariate continuous study - the one parameter grid
##########
# This file is where bugs B, C and D stop being possible.
#
# B  cts_miss_collect.R used mechanism = c("MAR", "AUX", "AUX-Y") while the DGM,
#    the analysis script and the check script all used "MNAR"/"MNAR-Y". Two of
#    the three mechanisms were collected from directories that never existed.
#    There is now one spelling: MNAR.
#
# C  cts_miss_metrics.R built its relative-efficiency reference from
#    method == "complete_data", but the collect grid omitted that method, so the
#    join found nothing and rel_efficiency was NA everywhere. "complete_data" is
#    in the grid below, so collect picks it up and the reference exists.
#
# D  cts_miss_analysis.R appended `filter(method == "complete_data")` AFTER
#    expand.grid, which renumbered every row. Array index i then meant one thing
#    to the analysis script and another to the check and rerun scripts, so a
#    failed_ids.txt would have resubmitted the wrong parameter rows. The grid is
#    now built once and never filtered; to run a subset, select on the VALUES
#    and keep the row numbers (see cts_miss_analysis.R).
#
# "complete_data" is the reference arm: no missingness is introduced at all, so
# the other methods can be scored against complete-data performance.

library(here)
source(here("R", "pipeline.R"))

MISS_METHODS_STUDY <- c("complete_cases", "mean_imputation", "missforest",
                        "regression", "missing_indicator", "IPW",
                        "multiple_imputation", "none", "complete_data")

grid <- expand.grid(
  scenario  = c(1, 2, 4, 5),
  n         = c(500),
  type      = c("both"),
  prop      = c(0.3),
  mechanism = c("MAR", "MNAR", "MNAR-Y"),
  method    = MISS_METHODS_STUDY,
  run       = c(1:100),
  stringsAsFactors = FALSE
)

# scenario 1 has no HTE, so a mechanism that depends on the treatment effect
# through U is not defined for it. Applied at construction, before any row
# numbers are handed out.
grid <- grid[!(grid$scenario == 1 & grid$mechanism == "MNAR-Y"), ]
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
