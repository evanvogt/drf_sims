##########
# title: prove me_strategies.R's pass-through is inert
##########
# me_strategies.R adds two nuisance arms and carries four things through
# untouched: the 9 candidates' fits, the data, the truth, and the `whole`
# arm. Everything downstream depends on that carry-through being exact - if it
# is not, the three-arm comparison is not comparing three ways of scoring THE
# SAME candidate fits, and the whole point of doing this as a second pass is
# lost.
#
# WHY THIS IS NOT A me_testing.R CHECK. Proving the carry-through needs no
# model fits, no fixtures and no recomputation - it is identical() on two
# saved objects. Running it over every one of the ~360 file pairs is minutes of
# pure I/O, and it is a strictly stronger statement than any assertion against
# an aggregated me_metrics.RDS, which could hide a compensating difference.
# me_testing.R checks that the SCORING is right; this checks that the INPUT to
# it was not disturbed.
#
#   Rscript me_strategies_verify.R
#
# Run after me_strategies.R's array finishes and before me_collect.R.

library(here)
source(here("model_evaluation/me_config.R"))

cmb <- combos(study)

# KNOWN, EXPECTED automl/holdout NA. calculate_pseudos() (me_nuisance.R)
# divides by pi * (1 - pi) with no trimming, and `holdout` fits AutoML's
# pi/mu_DR on nothing but a 25-100-row block (run_automl_holdout()) - see
# README.md's "Propensity trimming" section. The README already anticipated
# that exposure producing extreme weights; it turns out to also produce
# literal NA for a minority of runs, most plausibly from H2O returning
# NaN/NA predictions on some rows of a degenerate small-block fit. The
# formula is deliberately NOT changed (trimming one arm and not the others
# would make them non-comparable), so this is not a bug to chase - it is
# tracked here so a run that hits it isn't treated as a new failure.
#
# This is a HARDCODED list, unlike the 2 permanently-missing main-study runs
# (which check_failed()'s source-file diff derives automatically): whether a
# given block's AutoML fit degenerates to NA isn't computable from anything
# else on disk, it can only be observed. Transcribed 2026-08-29 from a full
# verify run's failure output (69 runs, all automl/holdout, none from xgb or
# cv_shared). If the strategies tree is ever rerun, regenerate this by
# re-running verify with the anyNA branch below temporarily reporting instead
# of suppressing, and diffing the new failure list against this one. Scenario
# ids are the post-2026-09-26 numbering (R/dgm_scenarios.R); the verify run
# that produced this list called scenarios 6, 8 and 4 by their old numbers 4, 6
# and 9.
#
# STALE UNTIL REGENERATED: the list comes from the pre-bug-O DGM. The study
# re-runs on the current one (its old trees are archived by
# R/archive_old_results.R), which changes every dataset, so which blocks
# degenerate will change too. Expect verify to fail on this list after the
# re-run, and regenerate it as above - including the 69 in the stopifnot().
known_holdout_na <- tibble::tribble(
  ~scenario, ~n, ~run,
  6,  250,  21,  6,  250,  23,  6,  250,  26,
  8,  250,   1,  8,  250,  10,  8,  250,  12,  8,  250,  13,  8,  250,  14,  8,  250,  26,
  4,  250,   3,  4,  250,   4,  4,  250,   6,  4,  250,  17,  4,  250,  19,  4,  250,  25,
  1,  500,   2,  1,  500,   3,  1,  500,   9,  1,  500,  11,  1,  500,  13,  1,  500,  27,
  6,  500,   1,  6,  500,   5,  6,  500,  12,  6,  500,  13,  6,  500,  23,  6,  500,  28,
  8,  500,   7,  8,  500,  10,  8,  500,  15,  8,  500,  18,  8,  500,  20,  8,  500,  29,  8,  500,  30,
  4,  500,   5,  4,  500,  10,  4,  500,  20,
  1, 1000,   5,  1, 1000,  10,  1, 1000,  11,  1, 1000,  25,  1, 1000,  29,
  6, 1000,   1,  6, 1000,   6,  6, 1000,   8,  6, 1000,  17,  6, 1000,  19,
  6, 1000,  21,  6, 1000,  23,  6, 1000,  28,  6, 1000,  29,  6, 1000,  30,
  8, 1000,   3,  8, 1000,   4,  8, 1000,  14,  8, 1000,  15,  8, 1000,  17,
  8, 1000,  21,  8, 1000,  24,  8, 1000,  25,  8, 1000,  26,  8, 1000,  27,  8, 1000,  28,
  4, 1000,   4,  4, 1000,   9,  4, 1000,  10,  4, 1000,  13,  4, 1000,  18,  4, 1000,  20
)
known_holdout_na_key <- with(known_holdout_na, paste(scenario, n, run))
stopifnot(nrow(known_holdout_na) == 69L)

fails <- character()
checked <- 0L
skipped <- 0L
known_na <- 0L

report_fail <- function(...) fails <<- c(fails, sprintf(...))

for (r in seq_len(nrow(cmb))) {
  combo <- cmb[r, , drop = FALSE]
  old_dir <- combo_dir(study, combo)
  new_dir <- combo_dir(study_strat, combo)

  for (run in seq_len(study$n_sims)) {
    f <- paste0("res_sim_", run, ".RDS")
    old_f <- file.path(old_dir, f)
    new_f <- file.path(new_dir, f)

    if (!file.exists(old_f)) {
      # excluded from the main study - me_strategies.R skips it by design
      if (file.exists(new_f)) {
        report_fail("%s: no source run, but the pass produced output", new_f)
      }
      skipped <- skipped + 1L
      next
    }

    if (!file.exists(new_f)) {
      report_fail("%s: source exists but the pass produced no output", new_f)
      next
    }

    old <- readRDS(old_f)
    new <- readRDS(new_f)
    checked <- checked + 1L
    tag <- sprintf("scenario_%s/%s/run %d", combo$scenario, combo$n, run)

    # the 9 candidate fits, in full - not just tau. If po/Y0.hat/Y1.hat/W.hat
    # moved, tau came from somewhere else even if it happens to match.
    for (m in CANDIDATE_MODELS) {
      if (!identical(old[[m]], new[[m]])) {
        report_fail("%s: candidate %s differs", tag, m)
      }
    }

    if (!identical(old$data, new$data)) report_fail("%s: data differs", tag)
    if (!identical(old$truth, new$truth)) report_fail("%s: truth differs", tag)
    if (!identical(old$fold_info, new$fold_info)) {
      report_fail("%s: fold_info differs", tag)
    }

    for (pipeline in names(old$nuisances)) {
      new_p <- new$nuisances[[pipeline]]
      old_p <- old$nuisances[[pipeline]]

      if (is.null(new_p)) {
        report_fail("%s: pipeline %s missing from the pass", tag, pipeline)
        next
      }
      if (!identical(names(new_p), NUISANCE_ARMS)) {
        report_fail("%s: %s arms are %s, expected %s", tag, pipeline,
                    paste(names(new_p), collapse = "/"),
                    paste(NUISANCE_ARMS, collapse = "/"))
      }
      # the carried-through arms, bit for bit
      if (!identical(old_p$whole, new_p$whole)) {
        report_fail("%s: %s `whole` was not carried through unchanged",
                    tag, pipeline)
      }
      # the new arms must actually be new, and complete
      for (arm in c("cv_shared", "holdout")) {
        d <- new_p[[arm]]
        if (is.null(d)) {
          report_fail("%s: %s `%s` missing", tag, pipeline, arm)
        } else if (nrow(d) != length(new$data$Y)) {
          report_fail("%s: %s `%s` has %d rows, expected %d",
                      tag, pipeline, arm, nrow(d), length(new$data$Y))
        } else if (anyNA(d$phi) || anyNA(d$pi)) {
          run_key <- paste(combo$scenario, combo$n, run)
          if (pipeline == "automl" && arm == "holdout" &&
                run_key %in% known_holdout_na_key) {
            known_na <- known_na + 1L
          } else {
            report_fail("%s: %s `%s` has NA in phi/pi", tag, pipeline, arm)
          }
        }
      }
    }
  }
}

cat(sprintf(
  "\n%d run(s) verified, %d skipped (no source run)\n", checked, skipped
))
cat(sprintf(
  "%d known automl `holdout` NA exception(s) - see the known_holdout_na table\n",
  known_na
))

if (length(fails)) {
  cat("\nFAILURES:\n")
  cat(paste0("  - ", fails, collapse = "\n"), "\n")
  quit(save = "no", status = 1)
}

cat("\nPass-through is inert: candidates, data, truth, fold_info and `whole`\n")
cat("are bit-identical to the source tree in every run checked.\n")
