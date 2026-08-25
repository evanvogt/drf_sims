##########
# title: prove me_strategies.R's pass-through is inert
##########
# me_strategies.R adds two nuisance arms and carries four things through
# untouched: the 9 candidates' fits, the data, the truth, and the `whole` and
# `cv_indep` (formerly `cv`) arms. Everything downstream depends on that
# carry-through being exact - if it is not, the four-arm comparison is not
# comparing four ways of scoring THE SAME candidate fits, and the whole point
# of doing this as a second pass rather than a rerun is lost.
#
# WHY THIS IS NOT A me_testing.R CHECK. Proving the carry-through needs no
# model fits, no fixtures and no recomputation - it is identical() on two
# saved objects. Running it over every one of the 358 file pairs is minutes of
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

fails <- character()
checked <- 0L
skipped <- 0L

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
      if (!identical(old_p$cv, new_p$cv_indep)) {
        report_fail("%s: %s `cv` -> `cv_indep` was not a pure rename",
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
        } else if (anyNA(d$phi) || anyNA(d$tau_T) || anyNA(d$pi)) {
          report_fail("%s: %s `%s` has NA in phi/tau_T/pi", tag, pipeline, arm)
        }
      }
      # cv_shared must use the CANDIDATES' folds, which is the entire
      # difference between it and cv_indep. If the two are identical, the
      # pass was handed the wrong fold vector.
      if (identical(new_p$cv_shared, new_p$cv_indep)) {
        report_fail("%s: %s `cv_shared` is identical to `cv_indep`",
                    tag, pipeline)
      }
    }
  }
}

cat(sprintf(
  "\n%d run(s) verified, %d skipped (no source run)\n", checked, skipped
))

if (length(fails)) {
  cat("\nFAILURES:\n")
  cat(paste0("  - ", fails, collapse = "\n"), "\n")
  quit(save = "no", status = 1)
}

cat("\nPass-through is inert: candidates, data, truth, fold_info, `whole` and\n")
cat("`cv_indep` are bit-identical to the source tree in every run checked.\n")
