##########
# title: correlated sample-size sets - summarise SuperLearner drops and failures
##########
# Every run saves the SuperLearner arm's dropped learners as results$sl_dropped
# (R/cate_models.R::cate_methods()): one row per (fold, model, learner), model
# being Y0 / Y1 (the per-arm outcome models), W (the propensity) or tau (stage
# 2), with the reason. Two kinds of row land there (R/sl_library.R):
#   pretest drops   pretest_superlearner() fits each candidate alone with 2-fold
#                   CV and drops one that errors or predicts non-finite values
#   "(whole fit)"   the live SuperLearner() errored, so every prediction in that
#                   fold/model is the training mean (mark_failed_fit()); or
#                   every ensemble weight was zero and the fit used the
#                   lowest-CV-risk learner alone (reason "all ensemble weights
#                   zero, ...", runs after 2026-10-09 only)
# Nothing downstream reads it. This summarises it for both outcomes. One more
# event is derived rather than saved: when the pretest drops every learner in a
# library, the fit falls back to SL.mean alone (bug K); that shows here as a
# fold/model whose drops cover the whole sl_libraries(n) library.
#
# NOT RECOVERABLE from the saved results - only in the PBS logs, as warnings:
# the all-zero-prediction fallback to the mean (sl_split_fit()), learners that
# pass the pretest but error inside the live 10-fold CV (weighted 0 by
# SuperLearner), and the kept learners' pretest warnings. Ensemble weights are
# not saved either.
#
# Reads the collected <prefix>_all.RDS, as the metrics scripts do. If it isn't
# there, reads the per-run files directly (get_results()), so this also runs
# before *_corr_collect.R.
#
# Writes <res_path>/<prefix>_sl_failures.RDS for each outcome, a list of:
#   runs       one row per run: folds, SL fits, and counts of pretest drops,
#              whole-fit failures and SL.mean fallbacks
#   drops      sl_dropped itself, one row per (run, fold, model, learner), with
#              reason_class
#   fallbacks  one row per (run, fold, model) whose pretest dropped every learner
#   cells      per (rho, scenario, n): runs and fits affected, by event
#   learners   per (rho, scenario, n, model, learner): fits dropped, runs affected
#   reasons    per (model, learner, reason): fits, runs, and where it happened
# and prints the cells, learners and reasons tables. Run from
# sample_size/correlated/:
#   Rscript corr_sl_failures.R                 # both outcomes
#   Rscript corr_sl_failures.R binary          # or just one

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
})
source(here("R", "metrics.R"))     # unnest_results()
source(here("R", "sl_library.R"))  # sl_libraries()

outcomes <- commandArgs(trailingOnly = TRUE)
if (length(outcomes) == 0) outcomes <- c("continuous", "binary")
stopifnot(all(outcomes %in% c("continuous", "binary")))

config_file <- c(continuous = "continuous/cts_corr_config.R",
                 binary     = "binary/bin_corr_config.R")

# sl_dropped's model column -> the sl_libraries() slot its library came from
LIB_SLOT <- c(Y0 = "Y", Y1 = "Y", W = "W", tau = "tau")
# fits per fold: Y0, Y1 and W in nuisance_sl(), tau in stage_2_sl()
FITS_PER_FOLD <- length(LIB_SLOT)

reason_class <- function(learner, reason) {
  case_when(
    learner == "(whole fit)" &
      startsWith(reason, "all ensemble weights zero") ~ "whole fit: zero weights",
    learner == "(whole fit)"               ~ "whole fit failed",
    startsWith(reason, "error:")           ~ "pretest: error",
    reason == "error inside SuperLearner"  ~ "pretest: error inside SuperLearner",
    reason == "non-finite predictions"     ~ "pretest: non-finite predictions",
    TRUE                                   ~ "other"
  )
}

#' One row per run, and its sl_dropped rows with the run's keys attached
#'
#' A run without dr_superlearner never ran the SuperLearner arm and is left out.
#' sl_dropped is NULL in a run where nothing was dropped.
sl_records <- function(study, all_results_df) {
  u <- unnest_results(study, all_results_df)
  recs <- lapply(seq_len(nrow(u$df)), function(i) {
    sim_res <- u$df$sim_res[[i]]
    if (is.null(sim_res$dr_superlearner)) return(NULL)
    key_row <- u$keys[i, , drop = FALSE]
    d <- sim_res$sl_dropped
    list(
      run = bind_cols(key_row, tibble(n_folds = length(unique(sim_res$fold_indices)))),
      drops = if (!is.null(d) && nrow(d) > 0) {
        bind_cols(key_row[rep(1, nrow(d)), , drop = FALSE], as_tibble(d))
      }
    )
  })
  recs <- recs[!vapply(recs, is.null, logical(1))]
  runs <- bind_rows(lapply(recs, `[[`, "run"))
  drops <- bind_rows(lapply(recs, `[[`, "drops"))
  if (nrow(drops) == 0) {
    # nothing dropped anywhere: keep the columns, so the summaries come out empty
    drops <- bind_cols(select(runs[0, ], -n_folds),
                       tibble(fold = integer(), model = character(),
                              learner = character(), reason = character()))
  }
  list(runs = runs,
       drops = mutate(drops, reason_class = reason_class(learner, reason)))
}

#' (run, fold, model) fits whose pretest dropped every learner in the library
#'
#' The library is sl_libraries(n) as it is now, so a run made under an older
#' library would be judged against the wrong one.
find_fallbacks <- function(drops, key_cols) {
  drops %>%
    filter(learner != "(whole fit)") %>%
    group_by(across(all_of(c(key_cols, "fold", "model")))) %>%
    summarise(dropped = list(unique(learner)), .groups = "drop") %>%
    filter(pmap_lgl(list(n, model, dropped), function(n, model, dropped) {
      all(sl_libraries(n)[[LIB_SLOT[[model]]]] %in% dropped)
    })) %>%
    select(-dropped)
}

summarise_sl <- function(study, all_results_df) {
  rec <- sl_records(study, all_results_df)
  path_cols <- study$path_cols
  key_cols <- c(path_cols, "run")
  drops <- rec$drops
  fallbacks <- find_fallbacks(drops, key_cols)

  # distinct fits with each event, so a fit with several dropped learners
  # counts once
  fits_with <- function(d, name) {
    d %>% distinct(across(all_of(c(key_cols, "fold", "model")))) %>%
      count(across(all_of(key_cols)), name = name)
  }

  runs <- rec$runs %>%
    mutate(n_fits = FITS_PER_FOLD * n_folds) %>%
    left_join(fits_with(filter(drops, learner != "(whole fit)"), "fits_pretest_drop"),
              by = key_cols) %>%
    left_join(fits_with(filter(drops, learner == "(whole fit)"), "fits_whole_fit"),
              by = key_cols) %>%
    left_join(fits_with(fallbacks, "fits_fallback"), by = key_cols) %>%
    mutate(across(c(fits_pretest_drop, fits_whole_fit, fits_fallback),
                  ~ coalesce(.x, 0L)))

  cells <- runs %>%
    group_by(across(all_of(path_cols))) %>%
    summarise(
      runs                  = n(),
      runs_pretest_drop     = sum(fits_pretest_drop > 0),
      runs_whole_fit        = sum(fits_whole_fit > 0),
      runs_fallback         = sum(fits_fallback > 0),
      fits                  = sum(n_fits),
      fits_pretest_drop     = sum(fits_pretest_drop),
      fits_whole_fit        = sum(fits_whole_fit),
      fits_fallback         = sum(fits_fallback),
      .groups = "drop"
    ) %>%
    mutate(pct_runs_pretest_drop = 100 * runs_pretest_drop / runs,
           pct_runs_whole_fit    = 100 * runs_whole_fit / runs)

  # a learner can be dropped at most once per (fold, model), so the share of
  # that model's fits in the cell is n_fits_dropped / sum(n_folds)
  folds_per_cell <- runs %>%
    group_by(across(all_of(path_cols))) %>%
    summarise(runs = n(), folds = sum(n_folds), .groups = "drop")
  learners <- drops %>%
    group_by(across(all_of(c(path_cols, "model", "learner")))) %>%
    summarise(fits_dropped = n(),
              runs_affected = n_distinct(run), .groups = "drop") %>%
    left_join(folds_per_cell, by = path_cols) %>%
    mutate(pct_fits = 100 * fits_dropped / folds,
           pct_runs = 100 * runs_affected / runs) %>%
    select(-folds) %>%
    arrange(across(all_of(path_cols)), model, desc(fits_dropped))

  reasons <- drops %>%
    group_by(model, learner, reason_class, reason) %>%
    summarise(fits = n(),
              runs = n_distinct(paste(!!!syms(key_cols))),
              where = paste(sort(unique(paste0(
                "rho ", rho, " sc", scenario, " n", n))), collapse = ", "),
              .groups = "drop") %>%
    arrange(desc(fits))

  list(runs = runs, drops = drops, fallbacks = fallbacks, cells = cells,
       learners = learners, reasons = reasons)
}

options(width = 200, dplyr.summarise.inform = FALSE)

for (outcome in outcomes) {
  source(here("sample_size", "correlated", config_file[[outcome]]))  # -> study

  collected <- file.path(study$res_path, paste0(study$prefix, "_all.RDS"))
  all_results_df <- if (file.exists(collected)) {
    readRDS(collected)
  } else {
    message(outcome, ": no ", basename(collected), "; reading the per-run files")
    get_results(study, workers = 1)
  }

  out <- summarise_sl(study, all_results_df)
  rm(all_results_df)
  gc()

  out_file <- file.path(study$res_path, paste0(study$prefix, "_sl_failures.RDS"))
  saveRDS(out, out_file)

  cat("\n========== ", outcome, " ==========\n", sep = "")
  cat(nrow(out$runs), " runs with a SuperLearner arm; ",
      sum(out$runs$fits_pretest_drop > 0), " with a pretest drop, ",
      sum(out$runs$fits_whole_fit > 0), " with a whole-fit failure, ",
      sum(out$runs$fits_fallback > 0), " with an SL.mean fallback\n", sep = "")

  cat("\n-- per cell --\n")
  print(as.data.frame(out$cells), digits = 3, row.names = FALSE)

  if (nrow(out$learners) > 0) {
    cat("\n-- dropped learners, per cell and model --\n")
    print(as.data.frame(out$learners), digits = 3, row.names = FALSE)

    cat("\n-- reasons --\n")
    print(as.data.frame(mutate(out$reasons, reason = substr(reason, 1, 90),
                               where = substr(where, 1, 60))),
          row.names = FALSE)
  }
  message(outcome, ": -> ", out_file)
}
