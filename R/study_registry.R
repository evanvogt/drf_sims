##########
# title: registry of every simulation study
##########
# One row per study: a display name and its config file, relative to the repo
# root. Sourced by check_all.R, tests/run_tests.R (jobscript ranges against
# each grid) and R/archive_old_results.R. To add a study, add a row; to retire
# one, delete it.

study_registry <- data.frame(
  study_name  = c(
    "competing_risk",
    "crossfitting",
    "crossfitting/confidence_intervals",
    "missing/continuous",
    "missing/binary",
    "missing/ci_example",
    "model_evaluation",
    "validation/continuous",
    "correlated/continuous",
    "correlated/binary",
    "correlated/confidence_intervals/continuous",
    "correlated/confidence_intervals/binary",
    "correlated/confidence_intervals/optimal_sf (cts)",
    "competing_risk/single_event",
    "validation/binary",
    "competing_risk/snr",
    "correlated/binary/snr"
  ),
  config_path = c(
    "competing_risk/surv_config.R",
    "crossfitting/cf_config.R",
    "crossfitting/confidence_intervals/cf_ci_config.R",
    "missing/continuous/cts_miss_config.R",
    "missing/binary/bin_miss_config.R",
    "missing/ci_example/cts_miss_ci_config.R",
    "model_evaluation/me_config.R",
    "validation/continuous/cts_val_config.R",
    "sample_size/continuous/cts_corr_config.R",
    "sample_size/binary/bin_corr_config.R",
    "sample_size/confidence_intervals/continuous/cts_corr_ci_config.R",
    "sample_size/confidence_intervals/binary/bin_corr_ci_config.R",
    "sample_size/confidence_intervals/optimal_sf/cts_corr_ci_sf_config.R",
    "competing_risk/single_event/se_config.R",
    "validation/binary/bin_val_config.R",
    "competing_risk/snr/snr_config.R",
    "sample_size/snr/bin_snr_config.R"
  ),
  stringsAsFactors = FALSE
)
