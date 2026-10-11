##########
# title: verification checks for the interim-analysis validation study, binary
##########
# Run this before submitting validation/binary/jobscripts/bin_val_1.sh:
#
#   Rscript validation/binary/bin_val_testing.R
#   Rscript validation/binary/bin_val_testing.R full   # adds check 8
#
# The helpers shared with continuous/ (TE-VIMs, TreeSHAP column order,
# interaction_pval's name-based indexing, the data-driven top-covariate tests)
# are checked by continuous/cts_val_testing.R; this covers what is
# binary-specific, and what a binary array needs before 3300 jobs go out.
#
# Checks, in order:
#   1. the packages a row needs load, including sandwich for the HC3 tests
#   2. the grid is 3300 rows - rho 0.5 only, scenarios 2-4, 11 interim
#      proportions, none of which stringifies to a float artefact
#      (as.character() is a directory name)
#   3. robust = TRUE gives HC3 p-values: they match a hand-computed sandwich
#      t-test, differ from the classical ones, are not small under a null, and
#      still return NA for a v with no contrast
#   4. the binary DGM, for each of scenarios 2-4: a trial of 1000 generates
#      (every risk inside [0, 1]), Y is 0/1, tau is a risk difference, and its
#      in-sample Var(tau) is printed. Scenario 2's tau is driven by X4 alone
#   5. split_trial(): one trial of 1000 splits n1 / 1000 - n1, and a run's chunk
#      1 at 0.25 is the first rows of its chunk 1 at 0.30 (nested, paired)
#   6. run_all_cate_methods on one 250-row binary chunk: all three estimators,
#      both measures, finite crossfit tau and TE-VIMs, and a step-by-step timing
#   7. chunk_validations_se() on check 6's fit: one row per standard-error kind,
#      each chunk_validations()'s own, and a binary Y's p_dd / p_proj plausible
#   8. (full only) one replicate end to end with bin_val_1.sh's "1 1"
#      arguments: the results path, the split, all four comparisons over all
#      three estimators and both standard-error kinds, plausible p-values

library(dplyr)
library(furrr)
library(grf)
library(here)

source(here("R", "utils.R"))
source(here("validation", "binary", "bin_val_dgms.R"))
source(here("validation", "binary", "bin_val_models.R"))
source(here("validation", "binary", "bin_val_config.R"))

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

models_expected <- c("causal_forest", "dr_random_forest", "dr_superlearner")

# =============================================================================
cat("\n=== 1. dependencies ===\n")

for (pkg in c("xgboost", "SHAPforxgboost", "sandwich",
              "SuperLearner", "glmnet", "gam", "earth", "ranger")) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  report(ok, sprintf("%s is installed%s", pkg,
                     if (ok) sprintf(" (%s)", as.character(packageVersion(pkg))) else ""))
}

# =============================================================================
cat("\n=== 2. the parameter grid ===\n")

props <- unique(study$grid$interim_prop)
report(nrow(study$grid) == 3300,
       sprintf("grid is 3300 rows (got %d)", nrow(study$grid)))
report(identical(sort(unique(study$grid$scenario)), 2:4) &&
         identical(unique(study$grid$rho), 0.5),
       "scenarios 2-4, rho = 0.5 only")
report(length(props) == 11 &&
         isTRUE(all.equal(sort(props), seq(0.25, 0.75, by = 0.05))),
       "11 interim proportions, 0.25 to 0.75 in 0.05 steps")
prop_strings <- as.character(sort(props))
report(all(nchar(prop_strings) <= 4),
       sprintf("no interim proportion stringifies to a float artefact (%s)",
               paste(prop_strings, collapse = ", ")))

# =============================================================================
cat("\n=== 3. HC3 interaction tests ===\n")

set.seed(201)
n_t <- 600
W_t <- rbinom(n_t, 1, 0.5)
v_t <- rbinom(n_t, 1, 0.5)
Y_t <- rbinom(n_t, 1, 0.3 + 0.05 * W_t + 0.1 * v_t + 0.15 * W_t * v_t)

fit_t <- lm(Y_t ~ W_t * v_t)
V_t <- sandwich::vcovHC(fit_t, type = "HC3")
t_hc3 <- coef(fit_t)[["W_t:v_t"]] / sqrt(V_t["W_t:v_t", "W_t:v_t"])
p_hand <- 2 * pt(-abs(t_hc3), df = fit_t$df.residual)
p_robust <- interaction_pval(Y_t, W_t, v_t, robust = TRUE)
p_classic <- interaction_pval(Y_t, W_t, v_t, robust = FALSE)

report(isTRUE(all.equal(p_robust, unname(p_hand))),
       sprintf("robust p-value is the hand-computed HC3 t-test (%.4g)", p_robust))
report(!isTRUE(all.equal(p_robust, p_classic)),
       sprintf("and differs from the classical one (%.4g)", p_classic))

set.seed(202)
Y_null <- rbinom(n_t, 1, 0.3 + 0.05 * W_t + 0.1 * v_t)
report(interaction_pval(Y_null, W_t, v_t, robust = TRUE) > 0.05,
       "a null interaction gives a non-significant robust p-value")
report(is.na(interaction_pval(Y_t, W_t, rep(0, n_t), robust = TRUE)),
       "a constant v still returns NA under robust = TRUE")

X_adj <- data.frame(a = rnorm(n_t), v = v_t)
report(is.finite(interaction_pval_adj(Y_t, W_t, X_adj, "v", robust = TRUE)),
       "interaction_pval_adj returns a finite robust p-value")

# =============================================================================
cat("\n=== 4. the binary DGM ===\n")

# generate_scenario_data() stops on any risk outside [0, 1], so a trial that
# generates at all has every risk inside
for (sc in sort(unique(study$grid$scenario))) {
  setup_rng_stream(1)
  gen_sc <- tryCatch(generate_binary_scenario_data(sc, 1000, study$grid$rho[1]),
                     error = function(e) { cat("  ", conditionMessage(e), "\n"); NULL })
  report(!is.null(gen_sc), sprintf("scenario %d: a trial of 1000 generates", sc))
  if (is.null(gen_sc)) next
  tau_sc <- gen_sc$truth$tau
  report(all(gen_sc$dataset$Y %in% c(0, 1)) && all(tau_sc > -1 & tau_sc < 1),
         sprintf("scenario %d: Y is 0/1 and tau a risk difference (range %.3f to %.3f)",
                 sc, min(tau_sc), max(tau_sc)))
  cat(sprintf("  NOTE  scenario %d: in-sample Var(tau) = %.5f (SD %.4f) at n = 1000\n",
              sc, var(tau_sc), sd(tau_sc)))
}

# scenario 2's trial is the one the split and one-chunk checks below use
setup_rng_stream(1)
gen_big <- generate_binary_scenario_data(2, 1000, study$grid$rho[1])
tau_t <- gen_big$truth$tau
report(isTRUE(all.equal(abs(cor(tau_t, tanh(gen_big$dataset$X4))), 1)),
       "scenario 2: tau is linear in tanh(X4) alone - X4 is the only modifier")

# =============================================================================
cat("\n=== 5. split_trial ===\n")

ch25 <- split_trial(gen_big, 0.25)
ch30 <- split_trial(gen_big, 0.30)
report(nrow(ch25$data1) == 250 && nrow(ch25$data2) == 750 &&
         nrow(ch25$truth1) == 250 && nrow(ch25$truth2) == 750,
       "0.25 splits 250 / 750, truth alongside")
report(isTRUE(all.equal(ch25$data1, ch30$data1[1:250, ])),
       "chunk 1 at 0.25 is the first 250 rows of chunk 1 at 0.30 (nested)")
report(identical(as.integer(rownames(ch25$data2)), seq_len(750)),
       "chunk 2's row names are reset")

# =============================================================================
cat("\n=== 6. run_all_cate_methods on a 250-row binary chunk ===\n")

sl_t <- system.time(
  fit_c <- run_all_cate_methods(data = ch25$data1, n_folds = chunk_folds(250),
                                verbose_timing = TRUE)
)[["elapsed"]]
cat(sprintf("  NOTE  one 250-row chunk took %.1f min with %d worker(s); by step:\n",
            sl_t / 60, workers))
print(round(unlist(fit_c$timings), 1))
fit_c$timings <- NULL

covars <- colnames(as.matrix(ch25$data1[, -c(1:2)]))
report(setequal(names(fit_c), models_expected),
       sprintf("all three estimators ran (got %s)", paste(names(fit_c), collapse = ", ")))
for (model in names(fit_c)) {
  m <- fit_c[[model]]
  report(length(m$tau) == 250 && all(is.finite(m$tau)),
         sprintf("%s: one finite tau per row", model))
  report(identical(colnames(m$te_vims), covars) && identical(colnames(m$shap_vims), covars) &&
           all(is.finite(unlist(m$te_vims[1, ]))) && all(is.finite(unlist(m$shap_vims[1, ]))),
         sprintf("%s: finite te_vims and shap_vims over the same covariates", model))
}

# =============================================================================
cat("\n=== 7. chunk_validations_se on a binary chunk ===\n")

# check 6's fit as both chunks - the bookkeeping, and that the new p-values are
# plausible with a 0/1 Y
d_se <- ch25$data1
v_se <- chunk_validations_se(fit_c, fit_c, d_se, d_se)
v_hc <- chunk_validations(fit_c, fit_c, d_se, d_se, robust = TRUE)
for (model in names(v_se$subgroups)) {
  sg <- v_se$subgroups[[model]]
  tv <- v_se$top_var_tests[[model]]
  report(nrow(sg) == 2 && setequal(sg$se, c("classical", "hc3")) &&
           isTRUE(all.equal(unname(unlist(sg[sg$se == "hc3", c("top", "bottom")])),
                            unname(v_hc$subgroups[[model]][c("top", "bottom")]))),
         sprintf("%s: subgroups, one row per kind, hc3 = chunk_validations(robust = TRUE)", model))
  p_new <- unlist(tv[c("p_dd", "p_proj")])
  report(nrow(tv) == 4 && all(is.na(p_new) | (p_new > 0 & p_new < 1)),
         sprintf("%s: top_var_tests has 4 rows and plausible p_dd / p_proj", model))
}

# =============================================================================
cat("\n=== 8. one replicate end to end ===\n")

if (!run_full) {
  cat("  SKIP  (pass 'full' to run - fits both chunks and writes a results file)\n")
} else {
  analysis <- here("validation", "binary", "bin_val_analysis.R")
  # "1 1" - one worker, one grf thread - is what bin_val_1.sh passes, so the
  # time reported below is the time an array job will take
  elapsed <- system.time(
    status <- system2("Rscript", c(shQuote(analysis), "1", "1", "1"),
                      stdout = NULL, stderr = NULL)
  )[["elapsed"]]

  report(status == 0, sprintf("bin_val_analysis.R 1 1 1 exits cleanly (%.1f min)",
                              elapsed / 60))
  cat(sprintf(paste0("  NOTE  a single replicate took %.1f min on one core; ",
                     "check against bin_val_1.sh's walltime\n"),
              elapsed / 60))

  # row 1 of the grid is interim_prop = 0.25
  param1 <- study$grid[1, ]
  out_file <- file.path(combo_dir(study, param1), "res_sim_1.RDS")
  report(file.exists(out_file), sprintf("results land at %s", out_file))

  if (file.exists(out_file)) {
    res <- readRDS(out_file)
    val <- res$validations

    n1 <- round(param1$n * param1$interim_prop)
    report(nrow(res$results1$data) == n1 && nrow(res$results2$data) == param1$n - n1,
           sprintf("the trial splits %d / %d", n1, param1$n - n1))
    report(setequal(names(val), c("subgroups", "variances", "var_imps", "top_var_tests")),
           "validations carries all four chunk comparisons")
    report(setequal(names(val$subgroups), models_expected),
           "every comparison covers all three estimators")

    plausible <- function(p) all(is.na(p) | (p > 0 & p < 1))
    p_cols <- c("p_cts", "p_cts_adj", "p_split", "p_dd", "p_proj")
    for (model in names(val$subgroups)) {
      sg <- val$subgroups[[model]]
      report(is.data.frame(sg) && setequal(sg$se, c("classical", "hc3")) &&
               plausible(c(sg$top, sg$bottom)),
             sprintf("%s: subgroup p-values, both kinds, are plausible (%s)", model,
                     paste(signif(c(sg$top, sg$bottom), 3), collapse = ", ")))
      tv <- val$top_var_tests[[model]]
      report(nrow(tv) == 4 && setequal(unique(tv$se), c("classical", "hc3")) &&
               all(p_cols %in% names(tv)) && plausible(unlist(tv[p_cols])),
             sprintf("%s: top_var_tests has both measures x both kinds and plausible p-values",
                     model))
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
cat("\nall checks passed. next: qsub validation/binary/jobscripts/bin_val_1.sh\n")
