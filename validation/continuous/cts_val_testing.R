##########
# title: verification checks for the interim-analysis validation study
##########
# Run this before submitting validation/continuous/jobscripts/cts_val_1.sh:
#
#   Rscript validation/continuous/cts_val_testing.R
#   Rscript validation/continuous/cts_val_testing.R full   # adds check 7
#
# Checks 1-6 are cheap and exercise the helpers directly. Check 7 runs one real
# replicate end to end and is worth the wait on the cluster, where it also tells
# you whether the jobscript's walltime is right before 1100 jobs go out.
#
# Checks, in order:
#   1. xgboost and SHAPforxgboost load. They are used nowhere else in this repo,
#      so a missing install here is the one thing that stops the array dead
#   2. the grid is 1100 rows over 11 interim proportions, and - the reason
#      cts_val_config.R rounds - none of those proportions stringifies to a
#      float artefact, since as.character(interim_prop) becomes a directory name
#   3. interaction_pval returns the W:v interaction, not the intercept. This is
#      the bottom_pval bug: the old code read coefficient row 1 instead of the
#      interaction, so the check is against both rows, not just the right one
#   4. interaction_pval returns NA on a covariate with no contrast, the case
#      where rpart predicted no bottom10 leaf into chunk 2 and the old
#      positional index would have read whatever row happened to be there
#   5. get_shap_vims returns one column per covariate IN X's column order.
#      shap.values() sorts its output by importance, and cts_val_analysis.R
#      ranks te_vims and shap_vims by position, so an unreindexed return would
#      scramble every rank silently rather than erroring
#   6. run_all_cate_methods attaches both te_vims and shap_vims to both
#      estimators, over the same covariates, so the rank comparison in
#      cts_val_analysis.R has matching rows to line up
#   7. (full only) one replicate end to end: results land in a directory named
#      by the un-mangled interim_prop, validations carries all four comparisons,
#      and bottom is a plausible p-value rather than the near-zero intercept
#      p-value the old code produced

library(dplyr)
library(furrr)
library(grf)
library(here)

source(here("R", "utils.R"))
source(here("validation", "continuous", "cts_val_dgms.R"))
source(here("validation", "continuous", "cts_val_models.R"))
source(here("validation", "continuous", "cts_val_config.R"))

args <- commandArgs(trailingOnly = TRUE)
run_full <- length(args) >= 1 && args[1] == "full"

workers <- 2

metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)

pass <- character()
fail <- character()

report <- function(ok, msg) {
  cat(if (ok) "  PASS  " else "  FAIL  ", msg, "\n", sep = "")
  if (ok) pass <<- c(pass, msg) else fail <<- c(fail, msg)
}

# =============================================================================
cat("\n=== 1. the new dependencies ===\n")

for (pkg in c("xgboost", "SHAPforxgboost")) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  report(ok, sprintf("%s is installed%s", pkg,
                     if (ok) sprintf(" (%s)", as.character(packageVersion(pkg))) else ""))
}

# =============================================================================
cat("\n=== 2. the parameter grid ===\n")

props <- unique(study$grid$interim_prop)

report(nrow(study$grid) == 1100,
       sprintf("grid is 1100 rows (got %d)", nrow(study$grid)))
report(length(props) == 11,
       sprintf("11 distinct interim proportions (got %d)", length(props)))
report(isTRUE(all.equal(sort(props), seq(0.25, 0.75, by = 0.05))),
       "interim proportions are 0.25 to 0.75 in 0.05 steps")

# as.character() is what file.path() uses to build the results directory name,
# so this is the check that cts_val_config.R's round() is doing its job
prop_strings <- as.character(sort(props))
report(all(nchar(prop_strings) <= 4),
       sprintf("no interim proportion stringifies to a float artefact (%s)",
               paste(prop_strings, collapse = ", ")))

# =============================================================================
cat("\n=== 3. interaction_pval reads the interaction, not the intercept ===\n")

# Built so the two p-values are nowhere near each other: a large intercept
# drives its p-value to ~0 whatever the interaction does, while the interaction
# effect is kept modest so its p-value stays well clear of the floor. Compared
# with identical() rather than all.equal() - both come from the same fit so
# exact equality holds, and all.equal's tolerance is absolute near zero, which
# would call a 1e-20 interaction p-value equal to an underflowed intercept one.
set.seed(101)
n_t <- 400
W_t <- rbinom(n_t, 1, 0.5)
v_t <- rbinom(n_t, 1, 0.5)
Y_t <- 50 + 0.5 * W_t + 0.5 * v_t + 0.45 * W_t * v_t + rnorm(n_t)

co_t <- summary(lm(Y_t ~ W_t * v_t))$coefficients
p_interaction <- unname(co_t[4, 4])
p_intercept <- unname(co_t[1, 4])
got <- interaction_pval(Y_t, W_t, v_t)

report(identical(got, p_interaction),
       sprintf("returns the W:v interaction p-value (%.3e)", p_interaction))
report(!identical(got, p_intercept),
       sprintf("does NOT return the intercept p-value (%.3e) - the old bug",
               p_intercept))
report(p_interaction > 1e-6 && p_intercept < p_interaction / 1e6,
       "the two p-values are far enough apart for that comparison to mean something")

# and with no interaction planted, the p-value must be free to be large - the
# intercept's never would be
set.seed(102)
Y_null <- 50 + 0.5 * W_t + 0.5 * v_t + rnorm(n_t)
report(interaction_pval(Y_null, W_t, v_t) > 0.05,
       "a null interaction gives a non-significant p-value")

# a continuous v takes the same path - this is what the top-variable test uses
x_cts <- rnorm(n_t)
Y_cts <- 50 + 0.5 * W_t + 2 * W_t * x_cts + rnorm(n_t)
report(interaction_pval(Y_cts, W_t, x_cts) < 0.01,
       "a planted continuous W x X interaction is detected")

# =============================================================================
cat("\n=== 4. interaction_pval's degenerate guard ===\n")

report(is.na(interaction_pval(Y_t, W_t, rep(0, n_t))),
       "constant v returns NA (rpart predicted no such leaf into chunk 2)")
report(is.na(interaction_pval(Y_t, W_t, rep(1, n_t))),
       "constant v = 1 returns NA")
report(!is.na(interaction_pval(Y_t, W_t, v_t)),
       "a v with contrast still returns a value")

# =============================================================================
cat("\n=== 5. get_shap_vims shape and column order ===\n")

set.seed(103)
n_s <- 300
X_s <- cbind(X1 = rnorm(n_s), X2 = rnorm(n_s), X4 = rnorm(n_s),
             X01 = rbinom(n_s, 1, 0.5))
tau_s <- -X_s[, "X4"] + rnorm(n_s, sd = 0.1)

sv <- get_shap_vims(X_s, tau_s)

report(nrow(sv) == 1 && ncol(sv) == ncol(X_s),
       sprintf("returns 1 x %d (got %d x %d)", ncol(X_s), nrow(sv), ncol(sv)))
report(identical(colnames(sv), colnames(X_s)),
       "columns are in X's order, not shap.values()' descending-importance order")
report(all(is.finite(unlist(sv[1, ]))) && all(unlist(sv[1, ]) >= 0),
       "all mean |SHAP| values are finite and non-negative")
report(colnames(sv)[which.max(unlist(sv[1, ]))] == "X4",
       "the planted CATE driver X4 is ranked most important")

# the order check above only bites if shap.values() really does return a
# different order - confirm the trap it guards against is live
raw_order <- names(SHAPforxgboost::shap.values(
  xgb_model = cvboost_cate(X_s, tau_s), X_train = X_s)$mean_shap_score)
report(!identical(raw_order, colnames(X_s)),
       sprintf("shap.values() really does reorder its output (%s) - the reindex is needed",
               paste(raw_order, collapse = ", ")))

# =============================================================================
cat("\n=== 6. run_all_cate_methods attaches both measures ===\n")

setup_rng_stream(1)
gen_t <- generate_continuous_scenario_data(3, 250)
fit_t <- run_all_cate_methods(data = gen_t$dataset, n_folds = 10)

covars <- colnames(as.matrix(gen_t$dataset[, -c(1:2)]))

for (model in names(fit_t)) {
  m <- fit_t[[model]]
  report(!is.null(m$te_vims) && !is.null(m$shap_vims),
         sprintf("%s: carries both te_vims and shap_vims", model))
  report(identical(colnames(m$te_vims), covars) &&
           identical(colnames(m$shap_vims), covars),
         sprintf("%s: both measures cover the same covariates in the same order", model))
  report(all(is.finite(unlist(m$shap_vims[1, ]))),
         sprintf("%s: shap_vims row is finite", model))
}

# =============================================================================
cat("\n=== 7. one replicate end to end ===\n")

if (!run_full) {
  cat("  SKIP  (pass 'full' to run - fits both chunks and writes a results file)\n")
} else {
  analysis <- here("validation", "continuous", "cts_val_analysis.R")
  elapsed <- system.time(
    status <- system2("Rscript", c(shQuote(analysis), "1"), stdout = NULL, stderr = NULL)
  )[["elapsed"]]

  report(status == 0, sprintf("cts_val_analysis.R 1 exits cleanly (%.1f min)",
                              elapsed / 60))
  cat(sprintf("  NOTE  a single replicate took %.1f min; jobscript walltime is 1h\n",
              elapsed / 60))

  # row 1 of the grid is interim_prop = 0.25
  out_file <- file.path(study$res_path, "scenario_3", "1000", "0.25", "res_sim_1.RDS")
  report(file.exists(out_file),
         sprintf("results land in a directory named by the un-mangled interim_prop (%s)",
                 out_file))

  if (file.exists(out_file)) {
    res <- readRDS(out_file)
    val <- res$validations

    report(setequal(names(val),
                    c("subgroups", "variances", "var_imps", "top_var_tests")),
           "validations carries all four chunk comparisons")

    for (model in names(val$subgroups)) {
      sg <- val$subgroups[[model]]
      # the old code put the intercept p-value here, which was essentially
      # always ~0 - a plausible p-value is the fix's end-to-end evidence
      report(is.na(sg[["bottom"]]) || (sg[["bottom"]] > 0 && sg[["bottom"]] < 1),
             sprintf("%s: bottom is a plausible p-value (%.3g), not the intercept's",
                     model, sg[["bottom"]]))
      report(is.na(sg[["top"]]) || (sg[["top"]] > 0 && sg[["top"]] < 1),
             sprintf("%s: top is a plausible p-value (%.3g)", model, sg[["top"]]))

      vi <- val$var_imps[[model]]
      report(nrow(vi) == 2 * length(covars) && "measure" %in% names(vi) &&
               setequal(unique(vi$measure), c("tevim", "shap")),
             sprintf("%s: var_imps has %d rows over both measures (got %d)",
                     model, 2 * length(covars), nrow(vi)))

      tv <- val$top_var_tests[[model]]
      report(nrow(tv) == 2 && setequal(unique(tv$measure), c("tevim", "shap")),
             sprintf("%s: top_var_tests has one row per measure", model))
      report(all(tv$x_top %in% covars) && all(tv$x_top2 %in% covars),
             sprintf("%s: x_top/x_top2 name real covariates (%s)",
                     model, paste(tv$x_top, collapse = ", ")))
      report(all(is.na(tv$p_cts) | is.finite(tv$p_cts)) &&
               all(is.na(tv$p_split) | is.finite(tv$p_split)),
             sprintf("%s: p_cts and p_split are finite or NA, never NaN", model))
    }
  }
}

# =============================================================================
plan(sequential)

cat("\n=== summary ===\n")
cat(sprintf("  %d passed, %d failed\n", length(pass), length(fail)))
if (length(fail) > 0) {
  cat("\nfailures:\n")
  for (f in fail) cat("  - ", f, "\n", sep = "")
  quit(status = 1)
}
cat("\nall checks passed. next: qsub validation/continuous/jobscripts/cts_val_1.sh\n")
