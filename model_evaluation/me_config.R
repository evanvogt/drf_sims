##########
# title: model evaluation study - the one definition of its parameter grid
##########
# Sourced by me_analysis.R, me_strategies.R, me_split.R, me_check.R,
# me_collect.R, me_metrics.R, me_strategies_verify.R and me_testing.R.
# The array index handed to me_analysis.R is a row number of `grid`, so `grid`
# must never be filtered or reordered after construction. The same applies to
# me_strategies.R, which is why study_strat below reuses this exact grid
# rather than one narrowed to the runs that completed.
#
# Scenario set is deliberately 4, not all 10: this study's research question
# (do cheap proxy losses rank 9 candidate models the way true PEHE would?)
# doesn't need every CATE-structure re-litigated, the same reasoning
# crossfitting/cf_analysis.R uses for its own scenario = c(1, 4, 6, 9).

library(here)
source(here("R", "pipeline.R"))

study <- study_config(
  name     = "model_evaluation",
  prefix   = "me",
  res_path = file.path(dirname(here()), "results", "model_evaluation"),
  grid = expand.grid(
    scenario = c(1, 4, 6, 9),
    n = c(250, 500, 1000),
    run = c(1:30),
    stringsAsFactors = FALSE
  ),
  path_cols   = c("scenario", "n"),
  n_sims      = 30,
  failed_file = here("model_evaluation", "jobscripts", "failed_ids.txt")
)

# The 9 candidate CATE-learner configurations, as returned by
# run_all_candidate_models() (me_models.R) and saved as top-level keys of
# each run's results (me_analysis.R). Declared once here so me_models.R and
# me_metrics.R cannot disagree about the model list - the same reason the
# grid above is declared once rather than re-typed in every script.
CANDIDATE_MODELS <- c("rf1", "rf2", "rf3", "net1", "net2", "net3", "SL1", "SL2", "SL3")

# ---- the nuisance-fitting arms ---------------------------------------------
# The four ways me_strategies.R estimates the evaluation nuisances, keyed under
# nuisances$<pipeline>$<arm> in each run's results. Declared once here for the
# same reason as CANDIDATE_MODELS: me_strategies.R writes these names and
# me_testing.R asserts against them, and a silent disagreement would show up
# only as the wrong number of score columns out of me_metrics.R.
#
# What varies across arms is what data the nuisance sees RELATIVE TO the data
# the candidate it is scoring was trained on. The 9 candidates are fit once,
# by me_analysis.R, single-crossfit over fold_info at V=10; every arm below
# scores those same fits.
#
#   whole      all n rows, fit and predicted           (no split at all)
#   cv_indep   90%, from a SECOND INDEPENDENT fold draw (retired - see below)
#   cv_shared  the candidate's own V-1 training folds, predicted on its held-out fold
#   holdout    the candidate's held-out fold only, fit AND predicted there
#
# cv_indep is retired but kept, because it is already computed and costs
# nothing to carry: me_analysis.R drew a second, independent fold assignment
# (nuisance_fold_info) on the theory that sharing one draw would correlate a
# candidate's tau_hat with the nuisance at the same row. That reasoning does
# not survive contact with the arithmetic. Row-level honesty - whether row i's
# own (Y_i, W_i) entered the model predicting at row i - holds under BOTH a
# shared and an independent draw. What the second draw changes is the OVERLAP
# of the two training sets, from 100% to (V-1)/V = 90%. It removes a tenth of
# the dependence, for a design that cannot be described from theory.
# crossfitting/README.md states the general form of this ("re-randomising the
# stage-2 split cannot remove that dependence"), which is why its scf_scf_new
# is an arm to be tested rather than the default. Keeping cv_indep alongside
# cv_shared makes that an empirical claim in this study too, rather than a
# deletion nobody can check.
NUISANCE_ARMS <- c("whole", "cv_indep", "cv_shared", "holdout")

# ---- the calibration score's group counts -----------------------------------
# How many tau_hat quantile groups the DR calibration score (calc_cal_score(),
# me_metrics.R) splits the evaluation rows into. Declared here for the same
# reason as NUISANCE_ARMS: me_metrics.R emits one column family per K and
# me_testing.R asserts the resulting column count, and a disagreement would
# surface only as the wrong number of score columns.
#
# BOTH values are carried rather than one being chosen, because the right K is
# a trade this study can measure instead of assume. At n = 250, K = 5 puts 50
# rows in each group - enough for a stable GATE^DR - while K = 10 puts 25,
# which resolves miscalibration more finely but estimates each group's mean
# pseudo-outcome from very little. Emitting both makes the sensitivity to K a
# result rather than a hidden choice.
CAL_QUANTILES <- c(5L, 10L)

# Blocks smaller than this are not worth handing to a 36-point XGBoost CV grid
# or a 20-model AutoML run: the `holdout` arm fits mu_DR/pi on nothing but the
# rows WITHIN one block, so a 25-row block leaves very little to fit on. At
# n = 250 the V = 10 candidate folds are 25 rows, so the arm pools adjacent
# fold PAIRS into 5 blocks of 50 (see holdout_blocks() in me_utils.R).
# Pooling rather than a fresh 5-fold draw is what keeps a block tied to the
# candidate split without refitting any candidate.
HOLDOUT_MIN_BLOCK <- 40L

# ---- the derived result trees ----------------------------------------------
# Both reuse `study`'s grid, path_cols and path_prefix, so combo_dir(),
# get_results(), check_failed() and compute_metrics() all work against them
# unchanged - the whole point of R/pipeline.R taking the study as an argument.

# me_strategies.R: a second pass over `study`'s existing results, adding the
# cv_shared and holdout arms. Same 360-row grid, because the array index must
# keep meaning the same row of the same grid. Written to a PARALLEL tree, not
# back into study$res_path - those 358 files are the only copy of a finished
# run.
study_strat <- study_config(
  name     = "model_evaluation_strategies",
  prefix   = "me_strat",
  res_path = file.path(dirname(here()), "results", "model_evaluation_strategies"),
  grid     = study$grid,
  path_cols   = study$path_cols,
  n_sims      = study$n_sims,
  failed_file = here("model_evaluation", "jobscripts", "failed_ids_strat.txt")
)

# me_split.R: the 80:20 arm. The only arm that refits the candidates - one fit
# on the 80% instead of a 10-fold crossfit - so it cannot share a results tree
# with anything above: its tau_hat exists only on the 20%, and its true_pehe
# has to be computed on those same rows to be a matched reference.
#
# n = 250 is excluded: a 50-row evaluation set is too thin to rank 9 models.
# The grid is filtered HERE, at construction, never afterwards - filtering a
# built grid renumbers every array index (bug D; see grid_indices()).
study_split <- study_config(
  name     = "model_evaluation_split",
  prefix   = "me_split",
  res_path = file.path(dirname(here()), "results", "model_evaluation_split"),
  grid = expand.grid(
    scenario = c(1, 4, 6, 9),
    n = c(500, 1000),
    run = c(1:30),
    stringsAsFactors = FALSE
  ),
  path_cols   = c("scenario", "n"),
  n_sims      = 30,
  failed_file = here("model_evaluation", "jobscripts", "failed_ids_split.txt")
)

# The fraction of each replicate the single-split arm trains the candidates on;
# the remainder is both the evaluation set and the only data its nuisance sees.
SPLIT_TRAIN_FRAC <- 0.8

#' The study object named by a CLI argument
#'
#' me_check.R, me_collect.R and me_metrics.R all serve three result trees now,
#' and each has to be told which. One lookup here rather than three copies of
#' the same switch().
#'
#' @param which "main" (the default), "strategies", or "split"
me_study <- function(which = "main") {
  switch(
    match.arg(which, c("main", "strategies", "split")),
    main       = study,
    strategies = study_strat,
    split      = study_split
  )
}
