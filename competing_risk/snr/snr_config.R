##########
# title: competing-risk SNR check - the one definition of its parameter grid
##########
# Sourced by snr_analysis.R, study.R (check, collect) and snr_metrics.R. The
# array index handed to snr_analysis.R is a row number of `grid`, so `grid`
# must never be filtered or reordered after construction.
#
# One factor at a time from a baseline (scenario 7, n = 500, hte_mult = 1),
# not factorial. Run varies fastest, so each cell is a contiguous block of 100:
#   rows   1-400  n = 500, hte_mult 1, 1.5, 2, 3  (jobscripts/snr_1.sh)
#   rows 401-600  n = 1000, 2000, hte_mult 1      (jobscripts/snr_2.sh)
#   rows 601-700  n = 4000, hte_mult 1            (jobscripts/snr_3.sh)
# Index 1 is the baseline run 1, 301 hte_mult 3 run 1, 601 n = 4000 run 1.
# The n blocks get their own jobscripts because rfsrc's memory grows roughly
# with n^2 (see README.md, Resources).
#
# A run is seeded by its run index alone (setup_rng_stream(run)), as in the
# main study, so the four n = 500 cells of a run share W, X and U and are
# paired. The n cells are not.
#
# Results live in results/competing_risk_snr/, not under results/competing_risk/:
# R/archive_old_results.R tars and removes that whole tree.

library(here)
source(here("R", "pipeline.R"))

# Fixed for every row: the main study's primary covariates (rho = 0.5) with
# censoring on, in its one scenario with HTE on both events
SNR_SCENARIO  <- 7
SNR_RHO       <- 0.5
SNR_CENSORING <- TRUE

# The production arms of the pseudo-value CF and the two DR-learners, plus the
# DR-learners' semi-oracle twins (true propensity 0.5). No SuperLearner, no
# scf/cvps arms, no ipw/csf_*: rel_eff_snr.R's comparison and nothing else.
SNR_MODELS <- c(
  "pseudo_cf_whole_oob",
  "pseudo_dr_whole_oob", "rsf_dr_oob",
  "pseudo_dr_whole_oob_semi", "rsf_dr_oob_semi"
)
# rel_eff_snr.R has an SNR for each event, none for the composite
SNR_ESTIMANDS <- c("RMTL1", "RMTL2")

SNR_RUNS <- 1:100

study <- study_config(
  name     = "competing_risk/snr",
  prefix   = "snr",
  res_path = file.path(dirname(here()), "results", "competing_risk_snr"),
  grid = rbind(
    expand.grid(run = SNR_RUNS, hte_mult = c(1, 1.5, 2, 3), n = 500),
    expand.grid(run = SNR_RUNS, n = c(1000, 2000, 4000), hte_mult = 1)
  )[, c("n", "hte_mult", "run")],
  path_cols   = c("n", "hte_mult"),
  path_prefix = c(hte_mult = "hte_"),
  n_sims      = length(SNR_RUNS),
  failed_file = here("competing_risk", "snr", "jobscripts", "failed_ids.txt")
)
