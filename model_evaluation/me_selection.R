##########
# title: model evaluation - per-run proxy selection and its references
##########
# The selection derivation shared by the thesis figures
# (results_processing/thesis_figures/model_eval.R) and the appendix tables
# (results_processing/thesis_tables/me_tables.R), so that a table cell and a
# figure's error bar come from the same runs. me_results.qmd keeps its own
# exploratory copy. Analysis only: display labels and layout stay with the
# calling script.
#
# Reads the strategies tree (me_metrics.R strategies): the DR and IF proxies,
# each with the estimated propensity or with it fixed at the true 0.5 (the
# `05` score types), crossed with the three fold regimes and the two nuisance
# pipelines. The calibration scores are left out.
#
# Source R/figures.R first: me_read_metrics() uses its label_factor(), and it
# loads dplyr and tidyr.

ME_KEY_COLS <- c("scenario", "n", "run", "model")
ME_N_CANDIDATES <- 9

# the score columns are <score_type>_<fold_type>_<pipeline>
# (model_evaluation/me_metrics.R)
ME_SCORE_TYPES <- c("dr", "dr05", "infl", "infl05")
ME_FOLD_TYPES <- c("cv_shared", "holdout", "whole")
ME_PIPELINES <- c("xgb", "automl")

ME_SCORE_COLS <- with(
  expand.grid(s = ME_SCORE_TYPES, f = ME_FOLD_TYPES, p = ME_PIPELINES,
              stringsAsFactors = FALSE),
  paste(s, f, p, sep = "_")
)

# anchored on the known tokens, so cv_shared's own underscore can't misparse
ME_SCORE_PATTERN <- paste0(
  "^(", paste(ME_SCORE_TYPES, collapse = "|"), ")_",
  "(", paste(ME_FOLD_TYPES, collapse = "|"), ")_",
  "(", paste(ME_PIPELINES, collapse = "|"), ")$"
)

#' Read the strategies tree's metrics, keeping the reported scenarios
#'
#' @param metrics_file path to me_strat_metrics.RDS
#' @param scenario_labels named by scenario number, in display order; other
#'   scenarios are dropped
#' @return one row per (scenario, n, run, model): true_pehe and the
#'   ME_SCORE_COLS, with scenario and n as factors in display order
me_read_metrics <- function(metrics_file, scenario_labels) {
  metrics <- readRDS(metrics_file)

  absent <- setdiff(ME_SCORE_COLS, names(metrics))
  if (length(absent) > 0) {
    stop(metrics_file, " has no ", paste(absent, collapse = ", "),
         ". Is it the strategies tree (Rscript me_metrics.R strategies)?",
         call. = FALSE)
  }

  metrics %>%
    filter(scenario %in% names(scenario_labels)) %>%
    select(all_of(ME_KEY_COLS), true_pehe, all_of(ME_SCORE_COLS)) %>%
    mutate(
      scenario = label_factor(scenario, scenario_labels),
      n = factor(n, levels = sort(unique(n)))
    )
}

#' One row per (scenario, n, run, proxy): how the proxy's pick did
#'
#' A candidate only enters a ranking if both its score and its truth are
#' there; regret is against the best candidate present in that run, as in
#' me_results.qmd.
#'
#' - hit, whether the pick is the candidate with the lowest true PEHE;
#' - pick_rank, the pick's true rank among the candidates present, 1 = best
#'   (ties take the lowest rank);
#' - regret, the pick's true PEHE minus the lowest true PEHE in the run.
#'
#' @param metrics output of me_read_metrics()
#' @return with the proxy's axes as id columns: score_type, fold_type,
#'   pipeline, and score_type split into score_family (dr, infl) and pi_regime
#'   (estimated, fixed)
me_selection <- function(metrics) {
  metrics %>%
    pivot_longer(all_of(ME_SCORE_COLS), names_to = "proxy_id",
                 values_to = "score") %>%
    filter(!is.na(score), !is.na(true_pehe)) %>%
    mutate(
      score_type = sub(ME_SCORE_PATTERN, "\\1", proxy_id),
      fold_type = sub(ME_SCORE_PATTERN, "\\2", proxy_id),
      pipeline = sub(ME_SCORE_PATTERN, "\\3", proxy_id)
    ) %>%
    group_by(scenario, n, run, score_type, fold_type, pipeline) %>%
    filter(n() >= 2) %>%
    summarise(
      hit = as.numeric(model[which.min(score)] == model[which.min(true_pehe)]),
      regret = true_pehe[which.min(score)] - min(true_pehe),
      pick_rank = rank(true_pehe, ties.method = "min")[which.min(score)],
      .groups = "drop"
    ) %>%
    mutate(
      score_family = sub("05$", "", score_type),
      pi_regime = ifelse(grepl("05$", score_type), "fixed", "estimated")
    )
}

#' Each candidate scored as if it were the pick, on me_selection()'s measures
#'
#' @param metrics output of me_read_metrics()
#' @return one row per (scenario, n, run, model) with a true PEHE: true_pehe,
#'   and the candidate's hit, pick_rank and regret in its run
me_candidates <- function(metrics) {
  metrics %>%
    filter(!is.na(true_pehe)) %>%
    group_by(scenario, n, run) %>%
    mutate(
      hit = as.numeric(model == model[which.min(true_pehe)]),
      pick_rank = rank(true_pehe, ties.method = "min"),
      regret = true_pehe - min(true_pehe)
    ) %>%
    ungroup() %>%
    select(all_of(ME_KEY_COLS), true_pehe, hit, pick_rank, regret)
}

#' The two references, per run, on the same measures as me_selection()
#'
#' Each candidate is scored as if it were the pick (me_candidates()). Then,
#' per (scenario, n):
#'
#' - random, a uniformly random pick: the run's mean over its candidates. Not
#'   a 0/1 indicator, so its hit takes the general MCSE;
#' - best_fixed, always deploying one candidate: the one with the lowest mean
#'   regret over runs. Chosen with hindsight from the truth, so an optimistic
#'   bound on any fixed rule. The same candidate is used for all three
#'   measures.
#'
#' @param metrics output of me_read_metrics()
#' @return one row per (scenario, n, run, ref), ref "random" or "best_fixed",
#'   with hit, pick_rank, regret, and `model` (the best fixed candidate; NA for
#'   random)
me_references <- function(metrics) {
  candidates <- me_candidates(metrics)

  random <- candidates %>%
    group_by(scenario, n, run) %>%
    summarise(across(c(hit, pick_rank, regret), mean), .groups = "drop") %>%
    mutate(ref = "random", model = NA_character_)

  best_model <- candidates %>%
    group_by(scenario, n, model) %>%
    summarise(value = mean(regret), .groups = "drop_last") %>%
    slice_min(value, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(scenario, n, model)

  best_fixed <- candidates %>%
    semi_join(best_model, by = c("scenario", "n", "model")) %>%
    transmute(scenario, n, run, hit, pick_rank, regret, ref = "best_fixed",
              model = as.character(model))

  bind_rows(random, best_fixed) %>%
    mutate(ref = factor(ref, levels = c("random", "best_fixed")))
}

#' Measures marginalised over the other design axes, per run
#'
#' Each level of each axis is averaged within run first, then (by the caller)
#' over runs, so its MCSE counts runs, not runs x the proxies sharing that
#' level (me_results.qmd pools the latter, which understates it). A hit
#' averaged this way is a proportion, so it takes the general MCSE.
#'
#' @param selection output of me_selection()
#' @param axes named list, one entry per axis in display order: `col`, the
#'   selection column it is read from, and `labels`, its levels' display
#'   labels named by value, in display order. Labels must be unique across axes
#' @param cols the measures to average
#' @return one row per (scenario, n, run, axis, level), axis and level factors
#'   in display order
me_axis_runs <- function(selection, axes, cols) {
  bind_rows(Map(function(axis, spec) {
    selection %>%
      mutate(axis = axis, level = unname(spec$labels[.data[[spec$col]]])) %>%
      group_by(scenario, n, run, axis, level) %>%
      summarise(across(all_of(cols), mean), .groups = "drop")
  }, names(axes), axes)) %>%
    mutate(
      axis = factor(axis, levels = names(axes)),
      level = factor(level,
                     levels = unlist(lapply(axes, function(a) unname(a$labels))))
    )
}
