##########
# title: correlated binary SNR check - the one definition of its parameter grid
##########
# The one grid definition: the array index is a row number of `grid`, so never
# filter or reorder it. Sourced by bin_snr_analysis.R, study.R and
# bin_snr_metrics.R.
#
# One factor at a time from n = 500, hte_mult = 1, run varying fastest:
#   rows    1-100   n = 500, scenario 1, hte_mult 1
#   rows  101-1300  n = 500, scenarios 2-4 x hte_mult 1, 1.5, 2, 3 (hte_mult faster)
#   rows 1301-2500  n = 1000, 2000, 4000 x scenarios 1-4 (scenario faster)
# Rows 1-1300 are jobscripts/bin_snr_1.sh, 1301-2500 bin_snr_2.sh.
#
# Every cell is the binary_snr set (binary_corr_0.5 with the treated risk
# clamped to [RD_EPS, 1 - RD_EPS]), bW calibrated at SNR_CALIB_N so the ATE is
# the same in every cell (README.md). Runs are seeded by run index, so the
# n = 500 cells share W and X; the n cells are not paired. Results go to
# results/correlated/binary_snr/.

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
