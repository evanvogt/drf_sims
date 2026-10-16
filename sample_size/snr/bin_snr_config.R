##########
# title: correlated binary SNR check - the one definition of its parameter grid
##########
# Sourced by bin_snr_analysis.R, study.R (check, collect) and
# bin_snr_metrics.R. The array index handed to bin_snr_analysis.R is a row
# number of `grid`, so `grid` must never be filtered or reordered after
# construction.
#
# One factor at a time from a baseline (n = 500, hte_mult = 1), not
# factorial, as competing_risk/snr/. Run varies fastest, so each cell is a
# contiguous block of 100:
#   rows    1-100   n = 500, scenario 1, hte_mult 1
#   rows  101-1300  n = 500, scenarios 2-4 x hte_mult 1, 1.5, 2, 3
#                   (hte_mult varies faster than scenario)
#   rows 1301-2500  n = 1000, 2000, 4000 x scenarios 1-4, hte_mult 1
#                   (scenario varies faster than n)
# Rows 1-1300 are jobscripts/bin_snr_1.sh, 1301-2500 bin_snr_2.sh. Index 101 is
# scenario 2 x1 run 1 (the main study's DGM), 1201 scenario 4 x3 run 1, 2101
# n = 4000 scenario 1 run 1.
#
# Every cell is binary_corr_0.5 (R/dgm_scenarios.R) as the "binary_snr" set,
# which clamps the treated risk to [RD_EPS, 1 - RD_EPS] once hte_mult pushes
# it out, with bW calibrated at SNR_CALIB_N whatever the n, so the ATE is the
# same in every cell (see README.md). A run is seeded by its run index alone,
# so the n = 500 cells of a run share W and X. The n cells are not paired.
#
# Results live in results/correlated/binary_snr/, beside the main study's
# results/correlated/binary/ rather than inside it.

library(here)
source(here("R", "pipeline.R"))

SNR_SET     <- "binary_snr"
SNR_CALIB_N <- 500

# the comparator and the random-forest DR-learner, with its semi-oracle (true
# propensity) and oracle (true outcome model too) versions. No SuperLearner.
# t_random_forest is derived from dr_random_forest's stage 1 at metrics time.
SNR_MODELS <- c("causal_forest", "dr_random_forest", "dr_oracle", "dr_semi_oracle")

SNR_RUNS <- 1:100

study <- study_config(
  name     = "correlated/binary/snr",
  prefix   = "bin_snr",
  res_path = file.path(dirname(here()), "results", "correlated", "binary_snr"),
  grid = rbind(
    expand.grid(run = SNR_RUNS, hte_mult = 1, scenario = 1, n = 500),
    expand.grid(run = SNR_RUNS, hte_mult = c(1, 1.5, 2, 3), scenario = 2:4, n = 500),
    expand.grid(run = SNR_RUNS, hte_mult = 1, scenario = 1:4, n = c(1000, 2000, 4000))
  )[, c("scenario", "n", "hte_mult", "run")],
  path_cols   = c("scenario", "n", "hte_mult"),
  path_prefix = c(scenario = "scenario_", hte_mult = "hte_"),
  n_sims      = length(SNR_RUNS),
  failed_file = here("sample_size", "snr", "jobscripts",
                     "failed_ids.txt")
)
