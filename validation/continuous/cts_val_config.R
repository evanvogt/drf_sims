##########
# title: interim-analysis validation study - the one definition of its parameter grid
##########
# Sourced by cts_val_analysis.R, study.R (check, collect) and
# cts_val_metrics.R. The array index handed to cts_val_analysis.R is a row
# number of `grid`, so `grid` must never be filtered or reordered after
# construction.
#
# Correlated covariates only, at rho = 0.5 - the correlated studies' primary arm
# (sample_size/, R/dgm_scenarios.R's CORRELATED COVARIATES). rho is a
# path column with the same "rho_" prefix those studies use, so results sit
# under rho_0.5/ and a rho = 0 arm would be a one-value change to the grid.
#
# Scenarios 1-4 (null, simple, complex, non-linear - R/dgm_scenarios.R,
# DESC_10): 4 x 11 x 100 = 4400 rows. In the null every significant interaction
# test is a type I error, so it is the size reference for the rest. Rows
# 4401-7100 take the null alone to 1000 runs at the three landmark interim
# points: 3 x 900 = 2700 rows. They are appended rather than merged in so rows
# 1-4400 keep their array indices; run is expand.grid's slowest column, so the
# first block is the grid exactly as it was before the extension.

library(here)
source(here("R", "pipeline.R"))

# round() is load-bearing, not cosmetic: interim_prop is a path_cols entry,
# so its as.character() form becomes a results directory name. Unrounded,
# seq(by = 0.05) yields values like 0.30000000000000004 and get_results()
# can no longer match the directory back to its grid row.
interim_props  <- round(seq(0.25, 0.75, by = 0.05), 2)
landmark_props <- c(0.25, 0.5, 0.75)

study <- study_config(
  name     = "validation/continuous",
  prefix   = "cts_val",
  res_path = file.path(dirname(here()), "results", "validation", "continuous"),
  grid = rbind(
    expand.grid(
      rho = 0.5,
      scenario = 1:4,
      n = 1000,
      interim_prop = interim_props,
      run = c(1:100),
      stringsAsFactors = FALSE
    ),
    # 1L, not 1: a double here would turn the whole scenario column double
    expand.grid(
      rho = 0.5,
      scenario = 1L,
      n = 1000,
      interim_prop = landmark_props,
      run = 101:1000,
      stringsAsFactors = FALSE
    )
  ),
  path_cols   = c("rho", "scenario", "n", "interim_prop"),
  path_prefix = c(rho = "rho_", scenario = "scenario_"),
  n_sims      = 100,
  failed_file = here("validation", "continuous", "jobscripts", "failed_ids.txt")
)
