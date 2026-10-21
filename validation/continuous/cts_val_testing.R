##########
# title: verification checks for the interim-analysis validation study
##########
# Run before submitting jobscripts/cts_val_1.sh:
#   Rscript validation/continuous/cts_val_testing.R
#   Rscript validation/continuous/cts_val_testing.R full   # adds check 9
# On the cluster, check 9's timing also sizes the jobscript's walltime.
#
#   1. xgboost, SHAPforxgboost and sandwich load (used nowhere else in the repo)
#   2. the 4400-row grid; no interim_prop stringifies to a float artefact (it
#      becomes a directory name); trial size and Var(tau)
#   3. interaction_pval returns the W:v row, not the intercept (the old
#      bottom_pval bug); interaction_pval_adj ignores a proxy of the modifier;
#      interaction_test()'s est/se/p and the tau_slope() contrasts
#   4. interaction_pval is NA for a covariate with no contrast
#   5. get_shap_vims keeps X's column order (shap.values() sorts by importance)
#   6. run_all_cate_methods fits all three estimators with both importance
#      measures over the same covariates; prints step timings
#   7. the data-driven top-covariate tests on a planted cos(x) modifier, a
#      binary covariate and a flat tau
#   8. chunk_validations_se(): classical and HC3 rows, the est_/se_/ref1_/true2_
#      columns and their sign conventions
#   9. (full) one replicate end to end with cts_val_1.sh's "1 1" arguments

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

# plus the DR SuperLearner's learners (R/sl_library.R). Those are in sim-env
# for the sample-size studies already, but a missing one would only show up as
# a learner the pretest silently drops
for (pkg in c("xgboost", "SHAPforxgboost", "sandwich",
              "SuperLearner", "glmnet", "gam", "earth", "ranger")) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  report(ok, sprintf("%s is installed%s", pkg,
                     if (ok) sprintf(" (%s)", as.character(packageVersion(pkg))) else ""))
}

# =============================================================================
cat("\n=== 2. the parameter grid ===\n")

props <- unique(study$grid$interim_prop)

report(nrow(study$grid) == 4400,
       sprintf("grid is 4400 rows (got %d)", nrow(study$grid)))
report(identical(sort(unique(study$grid$scenario)), 1:4) &&
         identical(unique(study$grid$rho), 0.5),
       "scenarios 1-4, rho = 0.5 only")
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

for (sc in sort(unique(study$grid$scenario))) {
  setup_rng_stream(1)
  gen_sc <- tryCatch(generate_continuous_scenario_data(sc, 1000, 0.5),
                     error = function(e) NULL)
  report(!is.null(gen_sc) && nrow(gen_sc$dataset) == 1000,
         sprintf("scenario %d generates a trial of 1000 at rho = 0.5", sc))
  if (!is.null(gen_sc)) {
    cat(sprintf("  NOTE  scenario %d: in-sample Var(tau) = %.4f\n", sc,
                var(gen_sc$truth$tau)))
  }
}

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

# the adjusted test: a proxy correlated with the real modifier interacts
# marginally but not once the modifier is in the model - the case p_cts_adj
# exists for under rho = 0.5
set.seed(104)
x_mod <- rnorm(n_t)
x_proxy <- 0.7 * x_mod + sqrt(1 - 0.7^2) * rnorm(n_t)
X_adj <- data.frame(Xmod = x_mod, Xproxy = x_proxy)
Y_adj <- 50 + 0.5 * W_t - 1 * W_t * x_mod + rnorm(n_t)
report(interaction_pval(Y_adj, W_t, x_proxy) < 0.01,
       "a proxy of the modifier interacts marginally (the trap)")
report(interaction_pval_adj(Y_adj, W_t, X_adj, "Xproxy") > 0.01,
       "interaction_pval_adj does not flag the proxy once the modifier is adjusted for")
report(interaction_pval_adj(Y_adj, W_t, X_adj, "Xmod") < 0.01,
       "interaction_pval_adj still detects the real modifier")

# the estimate and SE saved alongside the p-value: the classical
# triple is summary.lm's own W:v row, and the p-value is the one above
tr_t <- interaction_test(Y_t, W_t, v_t)
report(identical(unname(tr_t), unname(co_t[4, c(1, 2, 4)])) &&
         identical(names(tr_t), c("est", "se", "p")),
       sprintf("interaction_test returns summary.lm's W:v est/se/p (est %.3f)", tr_t[["est"]]))
report(identical(tr_t[["est"]], interaction_test(Y_t, W_t, v_t, robust = TRUE)[["est"]]),
       "the HC3 test has the same estimate, only its SE differs")
report(all(is.na(interaction_test(Y_t, W_t, rep(0, n_t)))),
       "a constant v gives an all-NA est/se/p")

# the reference contrasts: tau's projection slope on v, which for a 0/1 v is
# the difference in mean tau inside vs outside
tau_r <- 0.5 + 2 * v_t + rnorm(n_t, sd = 0.3)
report(isTRUE(all.equal(tau_slope(tau_r, v_t),
                        mean(tau_r[v_t == 1]) - mean(tau_r[v_t == 0]))),
       "tau_slope on a 0/1 v is the difference in mean tau")
report(is.na(tau_slope(NULL, v_t)) && is.na(tau_slope(tau_r, rep(1, n_t))),
       "tau_slope is NA without a tau or without a contrast")
tau_adj <- 0.5 - 1 * x_mod
report(isTRUE(all.equal(tau_slope_adj(tau_adj, X_adj, "Xmod"), -1)) &&
         abs(tau_slope_adj(tau_adj, X_adj, "Xproxy")) < 1e-8 &&
         tau_slope(tau_adj, x_proxy) < -0.5,
       "tau_slope_adj gives the proxy no slope once the modifier is in; tau_slope does")

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
gen_t <- generate_continuous_scenario_data(2, 250, study$grid$rho[1])
sl_t <- system.time(
  fit_t <- run_all_cate_methods(data = gen_t$dataset, n_folds = 5,
                                verbose_timing = TRUE)
)[["elapsed"]]
cat(sprintf("  NOTE  one 250-row chunk took %.1f min with %d worker(s); by step:\n",
            sl_t / 60, workers))
print(round(unlist(fit_t$timings), 1))
fit_t$timings <- NULL

covars <- colnames(as.matrix(gen_t$dataset[, -c(1:2)]))

report(setequal(names(fit_t), c("causal_forest", "dr_random_forest", "dr_superlearner")),
       sprintf("all three estimators ran (got %s)", paste(names(fit_t), collapse = ", ")))

sl_m <- fit_t$dr_superlearner
report(!is.null(sl_m) && length(sl_m$tau) == nrow(gen_t$dataset) && !anyNA(sl_m$tau),
       "dr_superlearner: one crossfit tau per row, none missing")
report(!is.null(sl_m) && all(is.finite(unlist(sl_m$te_vims[1, ]))),
       "dr_superlearner: every dropped-covariate refit filled its folds (finite TE-VIMs)")

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
cat("\n=== 7. the data-driven top-covariate tests ===\n")

# chunk 1 hands on cos(x) + noise; chunk 2's effect is cos(x). Symmetric in x,
# so a linear W x x term has nothing to find - scenario 4's situation
set.seed(105)
n_d1 <- 500
n_d2 <- 750
x_d1 <- rnorm(n_d1)
x_d2 <- rnorm(n_d2)
tau_d1 <- cos(x_d1) + rnorm(n_d1, sd = 0.2)
W_d2 <- rbinom(n_d2, 1, 0.5)
Y_d2 <- 0.5 * W_d2 + W_d2 * cos(x_d2) + rnorm(n_d2)

dd_t <- data_driven_split(x_d1, tau_d1, x_d2)
g_t <- projection_score(x_d1, tau_d1, x_d2)
report(length(dd_t$cuts) == 2 && dd_t$cuts[1] < 0 && dd_t$cuts[2] > 0,
       sprintf("a cos(x) CATE is split middle vs tails (cuts %s)",
               paste(signif(dd_t$cuts, 3), collapse = ", ")))
report(interaction_pval(Y_d2, W_d2, dd_t$high2, robust = TRUE) < 1e-3,
       "p_dd detects the cos(x) modifier")
report(!is.null(g_t) && interaction_pval(Y_d2, W_d2, g_t, robust = TRUE) < 1e-3,
       "p_proj detects the cos(x) modifier")
report(interaction_pval(Y_d2, W_d2, x_d2, robust = TRUE) > 0.05,
       "the linear p_cts does not - the scenario 4 blind spot")

# a binary covariate whose median is its upper value - X3 ~ Bernoulli(0.7)
xb_1 <- rbinom(n_d1, 1, 0.7)
xb_2 <- rbinom(n_d2, 1, 0.7)
report(all(as.numeric(xb_2 > median(xb_2)) == 0) &&
         identical(median_split(xb_2), as.numeric(xb_2)),
       "Bernoulli(0.7): the old median split is all 0, median_split() is the covariate")
dd_b <- data_driven_split(xb_1, 2 * xb_1 + rnorm(n_d1, sd = 0.2), xb_2)
report(identical(dd_b$high2, as.numeric(xb_2)),
       sprintf("and it is its own data-driven split (cut %s)",
               paste(dd_b$cuts, collapse = ", ")))
report(identical(projection_score(xb_1, tau_d1, xb_2), xb_2),
       "and its own projection score")

# a flat chunk 1 tau: no split, no curve, so NA rather than a test of nothing
dd_c <- data_driven_split(x_d1, rep(0.3, n_d1), x_d2)
report(length(dd_c$cuts) == 0 &&
         is.na(interaction_pval(Y_d2, W_d2, dd_c$high2, robust = TRUE)),
       "a flat chunk 1 tau gives no split and an NA p_dd")
report(is.null(projection_score(x_d1, rep(0.3, n_d1), x_d2)),
       "and no projection score, so an NA p_proj")

# =============================================================================
cat("\n=== 8. chunk_validations_se stacks both standard-error kinds ===\n")

# check 6's fit as both chunks - only the bookkeeping is under test here
d_se <- gen_t$dataset
v_se <- chunk_validations_se(fit_t, fit_t, d_se, d_se)
v_cl <- chunk_validations(fit_t, fit_t, d_se, d_se, robust = FALSE)
v_hc <- chunk_validations(fit_t, fit_t, d_se, d_se, robust = TRUE)
drop_se <- function(d) {
  d <- as.data.frame(d)
  d$se <- NULL
  rownames(d) <- NULL
  d
}
for (model in names(v_se$subgroups)) {
  sg <- v_se$subgroups[[model]]
  tv <- v_se$top_var_tests[[model]]
  same_sg <- all(vapply(list(list("classical", v_cl), list("hc3", v_hc)), function(k) {
    isTRUE(all.equal(unname(unlist(sg[sg$se == k[[1]], c("top", "bottom")])),
                     unname(k[[2]]$subgroups[[model]][c("top", "bottom")])))
  }, logical(1)))
  same_tv <- all(vapply(list(list("classical", v_cl), list("hc3", v_hc)), function(k) {
    isTRUE(all.equal(drop_se(tv[tv$se == k[[1]], ]),
                     drop_se(k[[2]]$top_var_tests[[model]])))
  }, logical(1)))
  report(nrow(sg) == 2 && same_sg && nrow(tv) == 4 && same_tv,
         sprintf("%s: classical and hc3 rows are chunk_validations()'s own", model))
}
report(identical(v_se$var_imps, v_hc$var_imps) &&
         identical(v_se$variances, v_cl$variances),
       "variances and var_imps do not depend on the standard errors")

# the estimates and reference contrasts. The truth is attached
# here as the analysis scripts attach it, so true2_* is filled
fit_tt <- fit_t
fit_tt$truth <- gen_t$truth
v_tt <- chunk_validations_se(fit_t, fit_tt, d_se, d_se)
forms <- c("cts", "cts_adj", "split", "dd", "proj")
sg_new <- paste0(rep(c("est", "se", "ref1", "true2"), each = 2), "_", c("top", "bottom"))
tv_new <- paste0(rep(c("est", "se", "ref1", "true2"), each = 5), "_", forms)
for (model in names(v_tt$subgroups)) {
  sg <- v_tt$subgroups[[model]]
  tv <- v_tt$top_var_tests[[model]]
  report(identical(names(sg), c("se", "top", "bottom", sg_new)) &&
           identical(names(tv)[-seq_len(11)], tv_new),
         sprintf("%s: est/se/ref1/true2 columns follow the p-values", model))
  same_kinds <- function(d, cols) {
    isTRUE(all.equal(d[d$se == "classical", cols], d[d$se == "hc3", cols],
                     check.attributes = FALSE))
  }
  report(same_kinds(sg, grep("^(est|ref1|true2)_", names(sg), value = TRUE)) &&
           same_kinds(tv, grep("^(est|ref1|true2)_", names(tv), value = TRUE)),
         sprintf("%s: est, ref1 and true2 are the same under both kinds", model))
  report(all(is.finite(unlist(sg[c("true2_top", "true2_bottom")]))) &&
           all(is.finite(unlist(tv[paste0("true2_", c("cts", "cts_adj", "split"))]))),
         sprintf("%s: true2 is filled once the truth is attached", model))
  # the sign conventions: chunk 1's top group has the larger estimated effect,
  # its bottom group the smaller, and the tree's high side is high
  report(all(sg$ref1_top >= 0 & sg$ref1_bottom <= 0, na.rm = TRUE) &&
           all(tv$ref1_dd >= 0, na.rm = TRUE),
         sprintf("%s: ref1_top >= 0, ref1_bottom <= 0, ref1_dd >= 0", model))
}
report(all(is.na(unlist(lapply(v_se$subgroups, function(d) d[c("true2_top", "true2_bottom")])))),
       "without the truth attached, true2 is NA")

# =============================================================================
cat("\n=== 9. one replicate end to end ===\n")

if (!run_full) {
  cat("  SKIP  (pass 'full' to run - fits both chunks and writes a results file)\n")
} else {
  analysis <- here("validation", "continuous", "cts_val_analysis.R")
  # "1 1" - one worker, one grf thread - is what cts_val_1.sh passes, so the
  # time reported below is the time an array job will take. Bare, the script
  # falls back to 5 workers and unthrottled grf and reports an optimistic time.
  elapsed <- system.time(
    status <- system2("Rscript", c(shQuote(analysis), "1", "1", "1"),
                      stdout = NULL, stderr = NULL)
  )[["elapsed"]]

  report(status == 0, sprintf("cts_val_analysis.R 1 1 1 exits cleanly (%.1f min)",
                              elapsed / 60))
  cat(sprintf(paste0("  NOTE  a single replicate took %.1f min on one core; ",
                     "check against cts_val_1.sh's walltime\n"),
              elapsed / 60))

  # row 1 of the grid is interim_prop = 0.25
  param1 <- study$grid[1, ]
  out_file <- file.path(combo_dir(study, param1), "res_sim_1.RDS")
  report(file.exists(out_file),
         sprintf("results land in a directory named by the un-mangled interim_prop (%s)",
                 out_file))

  if (file.exists(out_file)) {
    res <- readRDS(out_file)
    val <- res$validations

    n1 <- round(param1$n * param1$interim_prop)
    report(nrow(res$results1$data) == n1 &&
             nrow(res$results2$data) == param1$n - n1,
           sprintf("the trial splits %d / %d (got %d / %d)", n1, param1$n - n1,
                   nrow(res$results1$data), nrow(res$results2$data)))
    report(nrow(res$results1$truth) == nrow(res$results1$data) &&
             nrow(res$results2$truth) == nrow(res$results2$data),
           "each chunk's truth has one row per participant")

    report(setequal(names(val),
                    c("subgroups", "variances", "var_imps", "top_var_tests")),
           "validations carries all four chunk comparisons")
    report(setequal(names(val$subgroups),
                    c("causal_forest", "dr_random_forest", "dr_superlearner")),
           "every comparison covers all three estimators")

    plausible <- function(p) all(is.na(p) | (p > 0 & p < 1))
    for (model in names(val$subgroups)) {
      sg <- val$subgroups[[model]]
      report(is.data.frame(sg) && setequal(sg$se, c("classical", "hc3")),
             sprintf("%s: subgroups has one row per standard-error kind", model))
      # the old code put the intercept p-value here, which was essentially
      # always ~0 - a plausible p-value is the fix's end-to-end evidence
      report(plausible(sg$bottom),
             sprintf("%s: bottom is a plausible p-value (%s), not the intercept's",
                     model, paste(signif(sg$bottom, 3), collapse = ", ")))
      report(plausible(sg$top),
             sprintf("%s: top is a plausible p-value (%s)", model,
                     paste(signif(sg$top, 3), collapse = ", ")))

      vi <- val$var_imps[[model]]
      report(nrow(vi) == 2 * length(covars) && "measure" %in% names(vi) &&
               setequal(unique(vi$measure), c("tevim", "shap")),
             sprintf("%s: var_imps has %d rows over both measures (got %d)",
                     model, 2 * length(covars), nrow(vi)))

      tv <- val$top_var_tests[[model]]
      report(nrow(tv) == 4 && setequal(unique(tv$measure), c("tevim", "shap")) &&
               setequal(unique(tv$se), c("classical", "hc3")),
             sprintf("%s: top_var_tests has one row per measure x se", model))
      report(all(tv$x_top %in% covars) && all(tv$x_top2 %in% covars),
             sprintf("%s: x_top/x_top2 name real covariates (%s)",
                     model, paste(tv$x_top, collapse = ", ")))
      p_cols <- c("p_cts", "p_cts_adj", "p_split", "p_dd", "p_proj")
      report(all(c(p_cols, "dd_cuts", "dd_prop") %in% names(tv)) &&
               plausible(unlist(tv[p_cols])),
             sprintf("%s: %s are present and plausible or NA", model,
                     paste(p_cols, collapse = ", ")))

      # row 1 is scenario 1, the null: tau is constant, so every true contrast
      # is 0 (up to rounding) and every estimate is noise around it
      true2 <- unlist(c(sg[grep("^true2_", names(sg))], tv[grep("^true2_", names(tv))]))
      report(param1$scenario == 1 && length(true2) == 2 * 2 + 4 * 5 &&
               all(is.na(true2) | abs(true2) < 1e-8),
             sprintf("%s: scenario 1 (null) - every true2_* is 0", model))
      est <- unlist(c(sg[grep("^est_", names(sg))], tv[grep("^est_", names(tv))]))
      report(length(est) == 2 * 2 + 4 * 5 && all(is.na(est) | is.finite(est)),
             sprintf("%s: est_* present, finite or NA", model))
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
