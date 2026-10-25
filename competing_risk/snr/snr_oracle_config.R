##########
# title: competing-risk SNR check - the oracle DR-learner add-on
##########
# Sourced by snr_oracle_analysis.R, study.R (check, collect) and snr_metrics.R.
#
# The full-oracle DR-learner (true mu_0(x), mu_1(x) and e(x) = 0.5) fitted
# after the fact on the main SNR runs, so they need not be rerun: every
# res_sim_<run>.RDS of snr_config.R's study holds the run's data and truth, the
# pseudo-values are a deterministic function of the data (pseudo_all() draws no
# RNG), and the oracle needs no outcome model. Its own study, with the same
# grid, so row i here is row i of snr_config.R and the plumbing (check,
# rerun, collect) works unchanged. Results go to their own tree, so the main
# study's result files are never rewritten.
#
# One arm, not one per DR-learner: pseudo_dr_whole_oob and rsf_dr_oob use the
# same whole-sample pseudo-values in the correction term and differ only in
# mu-hat, so with the true mu they are the same estimator. The outcome in that
# term is still the pseudo-value, so the oracle removes nuisance error but
# keeps the pseudo-value noise and the stage-2 forest's error.

library(here)
source(here("competing_risk", "snr", "snr_config.R"))

# the main study, whose result files are read
snr_study <- study

ORACLE_ARM <- "pseudo_dr_oracle"

study <- study_config(
  name        = "competing_risk/snr/oracle",
  prefix      = "snr_oracle",
  res_path    = file.path(dirname(here()), "results", "competing_risk_snr_oracle"),
  grid        = snr_study$grid,
  path_cols   = snr_study$path_cols,
  path_prefix = snr_study$path_prefix,
  n_sims      = snr_study$n_sims,
  failed_file = here("competing_risk", "snr", "jobscripts", "failed_ids_oracle.txt")
)
