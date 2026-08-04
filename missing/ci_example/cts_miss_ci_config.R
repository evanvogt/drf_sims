##########
# title: MI confidence-interval example - the one parameter grid
##########
# A single cell of the missing-data design (MAR, multiple imputation), run over
# all five scenarios, asking how to build a CATE interval after imputation.
#
# The grid is arranged by (scenario, run) rather than left in expand.grid order.
# That was already true in cts_miss_ci_analysis.R and it is preserved here,
# because the PBS array index is a row number: reordering renumbers every job.

library(here)
source(here("R", "pipeline.R"))

grid <- expand.grid(
  scenario  = c(1:5),
  n         = c(500),
  type      = c("both"),
  prop      = c(0.3),
  mechanism = c("MAR"),
  method    = c("multiple_imputation"),
  run       = c(1:100),
  stringsAsFactors = FALSE
)
grid <- grid[order(grid$scenario, grid$run), ]
rownames(grid) <- NULL

study <- study_config(
  name     = "missing/ci_example",
  res_path = file.path(dirname(here()), "results", "missing", "ci_example"),
  grid        = grid,
  path_cols   = c("scenario", "n", "type", "prop", "mechanism", "method"),
  n_sims      = 100,
  failed_file = here("missing", "ci_example", "jobscripts", "failed_ids.txt")
)
