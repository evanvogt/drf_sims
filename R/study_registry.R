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
#                     scan (blocked = TRUE). No study is in this state right
#                     now - it's a mechanism kept for the next one that is.
#
# config_path is relative to the repo root (resolved with here::here()).
# config_var is the name study_config() is assigned to in that file - every
# current config uses "study", but the column exists so a future config that
# doesn't can still be added without changing check_all.R.
#
# patch_manifest names a one-off repair a study owes on top of being run, as the
# directory (relative to the study's res_path) that the repair writes its
# manifest into; NA for a study that owes none. It exists because "run" and
# "correct" are not the same state: missing/binary is 9,900/9,900 complete and
# was still owed the dr_random_forest HTE back-fill (R/patch_hte_tests.R), which
# the file counts alone report as "complete" and so hide. check_all.R reads the
# manifest rather than the result files - counting by opening 9,900 RDS objects
# would take half an hour a study and make the tracker too slow to run casually.
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
    "crossfit_rerun",
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
    "crossfitting strategy change - the last production study still double-crossfitting; now runs clean end-to-end, so this is its first run under the new strategy",
    "crossfitting strategy change",
    "crossfitting strategy change; also DGM bug A (continuous coefficients on logit scale)",
    "crossfitting strategy change",
    "crossfitting strategy change; also the DGM was wrong (see confidence_intervals/optimal_sf README)",
    "own comparison arms unchanged; only the production consumers of R/cate_models.R moved",
    "own comparison arms unchanged; pilot study, not part of the production rerun",
    "crossfitting strategy change; also bug F (dr_superlearner); plus the dr_random_forest HTE back-fill, patched in place - no re-run",
    "crossfitting strategy change; also the DGM was wrong three ways; plus the dr_random_forest HTE back-fill, patched in place - no re-run",
    "crossfitting strategy change",
    "first run, not a re-run - its own 9 candidates moved off double crossfitting (me_models.R); the 16 pre-change res_sim_*.RDS have been deleted, so the count restarts from zero",
    "crossfitting strategy change"
  ),
  # Only the two studies that ran under PROFILES$missing's old dr_rf_tests =
  # FALSE owe a patch. missing/ci_example is profile = "ci_mi" (tests off by
  # design) and has not started, so it will be born correct.
  patch_manifest = c(
    NA, NA, NA, NA, NA, NA, NA,
    NA, NA, "cts_miss_hte_patch", "bin_miss_hte_patch", NA, NA, NA
  ),
  blocked = c(
    FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
    FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
  ),
  stringsAsFactors = FALSE
)
