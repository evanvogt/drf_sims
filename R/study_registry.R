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
#   dgm_rerun       - needs resubmitting because the DGM changed under it
#                     (bug O, 2026-09-26); the results it had are archived
#                     by R/archive_old_results.R
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
    "sample_size/continuous/cts_config.R",
    "sample_size/binary/bin_config.R",
    "competing_risk/surv_config.R",
    "sample_size/confidence_intervals/continuous/cts_ci_config.R",
    "sample_size/confidence_intervals/binary/bin_ci_config.R",
    "sample_size/confidence_intervals/optimal_sf/cts_ci_sf_config.R",
    "sample_size/confidence_intervals/optimal_sf/bin_ci_sf_config.R",
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
    "dgm_rerun",
    "dgm_rerun",
    "crossfit_rerun",
    "crossfit_rerun",
    "crossfit_rerun",
    "dgm_rerun",
    "crossfit_rerun"
  ),
  reason = c(
    "crossfitting strategy change; also bug F (dr_superlearner)",
    "crossfitting strategy change; also bug F (dr_superlearner); also bug P and the risk-difference DGM (sample_size/binary/README.md)",
    "crossfitting strategy change - the last production study still double-crossfitting; now runs clean end-to-end, so this is its first run under the new strategy",
    "crossfitting strategy change",
    "crossfitting strategy change; also DGM bug A (continuous coefficients on logit scale); also bug P and the risk-difference DGM (sample_size/binary/README.md)",
    "crossfitting strategy change",
    "crossfitting strategy change; also the DGM was wrong (see sample_size/confidence_intervals/optimal_sf README); also bug P and the risk-difference DGM (sample_size/binary/README.md)",
    "own comparison arms unchanged by the crossfitting change, but bug O changed the continuous DGM it runs on, and the 2026-09-26 renumbering its scenario ids (1/4/6/9 -> 1/4/6/8); old results archived",
    "pilot study, not part of the production rerun; re-run because bug O changed the continuous DGM and the 2026-09-26 renumbering its scenario ids (1/6/9 -> 1/4/8); old results archived",
    "crossfitting strategy change; also bug F (dr_superlearner); plus bug O (continuous DGM) - all 12,600 rows re-run; plus scenario 2 (added later; was numbered 6 before the 2026-09-26 renumbering), rows 9901-12600",
    "crossfitting strategy change; also the DGM was wrong three ways; plus scenario 2 (added later; was numbered 6 before the 2026-09-26 renumbering), rows 9901-12600; plus bug N (MNAR-Y truth), repaired at metrics time - no re-run; plus bug P and the risk-difference DGM - all 12,600 rows re-run, which also makes bug N moot",
    "crossfitting strategy change",
    "its 358/360 runs (the first under single crossfitting, me_models.R) predate bug O, which changed the continuous DGM, and the 2026-09-26 renumbering (1/4/6/9 -> 1/4/6/8); archived with the strategies and split trees, so the count restarts from zero",
    "crossfitting strategy change"
  ),
  blocked = c(
    FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
    FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
  ),
  stringsAsFactors = FALSE
)
