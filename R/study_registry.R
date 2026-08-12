##########
# title: registry of every simulation study, for the unified rerun checker
##########
# One row per study. Sourced only by check_all.R (repo root) - none of the
# per-study *_check.R / *_collect.R / *_metrics.R scripts need this, so
# adding or removing a study here can never affect the existing per-study
# rerun workflow.
#
# category is one of:
#   crossfit_rerun - needs resubmitting because of the crossfitting strategy
#                     change, possibly plus a study-specific bug fix (see
#                     `reason`)
#   first_run       - never successfully run before; not a rerun
#   no_rerun        - already fine as-is; own arms/results are unaffected
#   broken          - currently fails to run; excluded from the filesystem
#                     scan (blocked = TRUE)
#
# config_path is relative to the repo root (resolved with here::here()).
# config_var is the name study_config() is assigned to in that file - every
# current config uses "study", but the column exists so a future config that
# doesn't can still be added without changing check_all.R.
#
# To add a new study: add one row. To retire a study: delete its row.

study_registry <- data.frame(
  study_name  = c(
    "continuous",
    "binary",
    "competing_risk",
    "confidence_intervals/continuous",
    "confidence_intervals/binary",
    "confidence_intervals/optimal_sf (cts)",
    "confidence_intervals/optimal_sf (bin)",
    "crossfitting",
    "crossfitting/confidence_intervals",
    "missing/continuous",
    "missing/binary",
    "missing/ci_example",
    "model_evaluation",
    "validation/continuous"
  ),
  config_path = c(
    "continuous/cts_config.R",
    "binary/bin_config.R",
    "competing_risk/surv_config.R",
    "confidence_intervals/continuous/cts_ci_config.R",
    "confidence_intervals/binary/bin_ci_config.R",
    "confidence_intervals/optimal_sf/cts_ci_sf_config.R",
    "confidence_intervals/optimal_sf/bin_ci_sf_config.R",
    "crossfitting/cf_config.R",
    "crossfitting/confidence_intervals/cf_ci_config.R",
    "missing/continuous/cts_miss_config.R",
    "missing/binary/bin_miss_config.R",
    "missing/ci_example/cts_miss_ci_config.R",
    "model_evaluation/me_config.R",
    "validation/continuous/cts_val_config.R"
  ),
  config_var = "study",
  category = c(
    "crossfit_rerun",
    "crossfit_rerun",
    "broken",
    "crossfit_rerun",
    "crossfit_rerun",
    "crossfit_rerun",
    "crossfit_rerun",
    "no_rerun",
    "no_rerun",
    "crossfit_rerun",
    "crossfit_rerun",
    "crossfit_rerun",
    "first_run",
    "crossfit_rerun"
  ),
  reason = c(
    "crossfitting strategy change; also bug F (dr_superlearner)",
    "crossfitting strategy change; also bug F (dr_superlearner)",
    "currently fails to run - independent estimator code, untouched by the crossfitting strategy change",
    "crossfitting strategy change",
    "crossfitting strategy change; also DGM bug A (continuous coefficients on logit scale)",
    "crossfitting strategy change",
    "crossfitting strategy change; also the DGM was wrong (see confidence_intervals/optimal_sf README)",
    "own comparison arms unchanged; only the production consumers of R/cate_models.R moved",
    "own comparison arms unchanged; pilot study, not part of the production rerun",
    "crossfitting strategy change; also bug F (dr_superlearner)",
    "crossfitting strategy change; also the DGM was wrong three ways",
    "crossfitting strategy change",
    "first run, not a re-run - but its own 9 candidates moved off double crossfitting (me_models.R), so its first 16 res_sim_*.RDS predate that and need deleting before continuing",
    "crossfitting strategy change"
  ),
  blocked = c(
    FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE,
    FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
  ),
  stringsAsFactors = FALSE
)
