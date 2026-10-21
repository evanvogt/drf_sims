##########
# title: prove me_strategies.R's pass-through is inert
##########
# me_strategies.R fits all three nuisance arms and carries four things
# through untouched: the 9 candidates' fits, the data, the truth, and
# fold_info. Everything downstream depends on that carry-through being exact -
# if it is not, the three-arm comparison is not comparing three ways of
# scoring the same candidate fits, and the whole point of doing this as a
# second pass is lost. `whole` is recomputed rather than carried
# (me_strategies.R's header), so it is checked for completeness like the other
# two arms, not for identity with the main tree.
#
# WHY THIS IS NOT A me_testing.R CHECK. Proving the carry-through needs no
# model fits, no fixtures and no recomputation - it is identical() on two
# saved objects. Running it over every one of the ~360 file pairs is minutes of
# pure I/O, and it is a strictly stronger statement than any assertion against
# an aggregated me_metrics.RDS, which could hide a compensating difference.
# me_testing.R checks that the scoring is right; this checks that the input to
# it was not disturbed.
#
#   Rscript me_strategies_verify.R
#
# Run after me_strategies.R's array finishes and before me_collect.R.

library(here)
source(here("R", "utils.R"))                      # trim_ps()
source(here("model_evaluation/me_config.R"))
source(here("model_evaluation", "me_nuisance.R")) # calculate_pseudos()

cmb <- combos(study)

# WHAT "COMPLETE" MEANS FOR AN ARM. The saved phi is not what gets scored:
# me_metrics.R recomputes it from mu0_DR/mu1_DR/pi with calculate_pseudos(),
# which trims pi to [0.05, 0.95] (DEVLOG bug S). So an arm is complete when
# its predictions have no NA and the pseudo-outcomes recomputed from them are
# all finite - checked here the same way the scoring computes them. A file
# written before the trimming can hold a NaN/Inf saved phi and still pass.
#
# Propensities of exactly 0 or 1 are not failures (the trimming is what
# handles them) but they are tallied per pipeline and arm, as the record of
# how often the trimming is what kept a run in. Before it, AutoML did this in
# about a quarter of holdout runs, 5% of whole and 0.5% of cv_shared, and
# XGBoost never; that replaced a hard-coded list of excused automl/holdout
# runs, which had gone stale with the DGM and missed the other two arms.

fails <- character()
checked <- 0L
skipped <- 0L
# runs with at least one pi of exactly 0 or 1, per pipeline/arm
pi_boundary <- integer()

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

    for (pipeline in c("xgb", "automl")) {
      new_p <- new$nuisances[[pipeline]]

      if (is.null(new_p)) {
        report_fail("%s: pipeline %s missing from the pass", tag, pipeline)
        next
      }
      if (!identical(names(new_p), NUISANCE_ARMS)) {
        report_fail("%s: %s arms are %s, expected %s", tag, pipeline,
                    paste(names(new_p), collapse = "/"),
                    paste(NUISANCE_ARMS, collapse = "/"))
      }
      # every arm is fitted by the pass, so every arm must be complete
      for (arm in NUISANCE_ARMS) {
        d <- new_p[[arm]]
        if (is.null(d)) {
          report_fail("%s: %s `%s` missing", tag, pipeline, arm)
        } else if (nrow(d) != length(new$data$Y)) {
          report_fail("%s: %s `%s` has %d rows, expected %d",
                      tag, pipeline, arm, nrow(d), length(new$data$Y))
        } else if (anyNA(d[c("mu0_DR", "mu1_DR", "pi")])) {
          # trimming cannot repair a missing prediction
          report_fail("%s: %s `%s` has NA in mu0_DR/mu1_DR/pi",
                      tag, pipeline, arm)
        } else {
          scored <- calculate_pseudos(d, new$data$Y, new$data$W)
          if (!all(is.finite(scored$phi))) {
            report_fail("%s: %s `%s` has non-finite phi even with pi trimmed",
                        tag, pipeline, arm)
          }
          if (any(d$pi %in% c(0, 1))) {
            key <- paste(pipeline, arm)
            pi_boundary[key] <- sum(pi_boundary[key], 1L, na.rm = TRUE)
          }
        }
      }
    }
  }
}

cat(sprintf(
  "\n%d run(s) verified, %d skipped (no source run)\n", checked, skipped
))
cat("\nruns with a propensity of exactly 0 or 1, kept in by the trimming:\n")
if (length(pi_boundary)) {
  for (key in sort(names(pi_boundary))) {
    cat(sprintf("  %-20s %d\n", key, pi_boundary[[key]]))
  }
} else {
  cat("  none\n")
}

if (length(fails)) {
  cat("\nFAILURES:\n")
  cat(paste0("  - ", fails, collapse = "\n"), "\n")
  quit(save = "no", status = 1)
}

cat("\nPass-through is inert: candidates, data, truth and fold_info are\n")
cat("bit-identical to the source tree, and all three arms are complete (no NA\n")
cat("prediction, finite trimmed phi), in every run checked.\n")
