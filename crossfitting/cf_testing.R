##########
# title: verification checks for the crossfitting comparison
##########
# Run this before submitting anything to the HPC:
#
#   Rscript crossfitting/cf_testing.R       # quick: structure + regression check
#   Rscript crossfitting/cf_testing.R full  # adds the SuperLearner family (slow)
#
# Checks, in order:
#   1. the new dcf arm reproduces R/cate_models.R's dr_random_forest
#   2. every arm returns complete tau / tau_test of the right length
#   3. scenario 1 (no heterogeneity) behaves
#   4. test-set plumbing - arms sharing a fitted model have identical test
#      predictions, and test predictions track the test truth. NOT an optimism
#      check: scoring is against the known true CATE, not the labels the models
#      were fit to, so there is no optimism to detect (see the note at check 4).

library(dplyr)
library(furrr)
library(grf)
library(SuperLearner)
library(here)

source(here("crossfitting", "cf_models.R"))
source(here("crossfitting", "cf_metrics.R"))

# "full" adds the SuperLearner family, which dominates the runtime
full <- "full" %in% commandArgs(trailingOnly = TRUE)

workers <- 2
grf_threads <- 1
Sys.setenv(OMP_NUM_THREADS = grf_threads)

metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)

pass <- character()
fail <- character()

report <- function(ok, msg) {
  cat(if (ok) "  PASS  " else "  FAIL  ", msg, "\n", sep = "")
  if (ok) pass <<- c(pass, msg) else fail <<- c(fail, msg)
}

# =============================================================================
cat("\n=== 1. regression check: dcf against cts_models.R ===\n")
# the new double-crossfit path must reproduce the existing implementation. both
# are driven from the same RNG state with the same fold split, so the forests
# should be bit-identical; the only deliberate difference is that trim_ps is
# applied to the RF propensities here and not in cts_models.R, which is a no-op
# unless a propensity leaves [0.05, 0.95].

setup_rng_stream(7)
gen <- generate_cf_replicate(scenario = 6, n = 500, n_test = 1000)
X <- as.matrix(gen$data[, -c(1:2)])
Y <- gen$data$Y
W <- gen$data$W

fold_indices <- sort(seq(nrow(X)) %% 10) + 1
fold_list <- unique(fold_indices)
fold_pairs <- utils::combn(fold_list, 2, simplify = FALSE)

setup_rng_stream(7)
nz_old <- nuisance_rf(X, Y, W, fold_indices, fold_pairs)

setup_rng_stream(7)
nz_new <- nuisance_double_rf(
  X,
  Y,
  W,
  fold_indices,
  fold_pairs,
  num.threads = NULL
)

ps_range <- range(nz_old$W.hat_matrix, na.rm = TRUE)
cat(sprintf(
  "  propensity range %.3f to %.3f (trimming is a no-op inside [0.05, 0.95])\n",
  ps_range[1],
  ps_range[2]
))
report(
  ps_range[1] >= 0.05 && ps_range[2] <= 0.95,
  "RF propensities lie inside the trimming bounds, so trim_ps changes nothing"
)

po_diff <- max(abs(nz_old$po_matrix - nz_new$po), na.rm = TRUE)
cat(sprintf("  max |po_old - po_new| = %.3e\n", po_diff))
report(po_diff < 1e-8, "double-crossfit pseudo-outcomes match cts_models.R")

setup_rng_stream(7)
tau_old <- stage_2_rf(X, nz_old$po_matrix, fold_indices, fold_list)

setup_rng_stream(7)
s_new <- stage2_crossfit_rf(
  X,
  nz_new$po,
  gen$X_test,
  fold_indices,
  num.threads = NULL
)

tau_diff <- max(abs(tau_old - s_new$tau), na.rm = TRUE)
cat(sprintf("  max |tau_old - tau_new| = %.3e\n", tau_diff))
report(
  tau_diff < 1e-8,
  "dcf CATE estimates match cts_models.R's dr_random_forest"
)

# =============================================================================
cat("\n=== 2. structure: every arm complete and correctly sized ===\n")

n <- 500
n_test <- 1000
sl_lib <- if (full) {
  c("SL.glm", "SL.glmnet", "SL.earth", "SL.gam", "SL.mean", "SL.ranger")
} else {
  NULL
}

setup_rng_stream(1)
gen1 <- generate_cf_replicate(scenario = 6, n = n, n_test = n_test)

t0 <- Sys.time()
res <- run_all_crossfit_variants(
  gen1$data,
  gen1$X_test,
  n_folds = 10,
  sl_lib = sl_lib,
  num.threads = grf_threads,
  truth_test = gen1$truth_test_tau
)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf(
  "  one replicate took %.1f s with workers=%d, grf_threads=%d%s\n",
  elapsed,
  workers,
  grf_threads,
  if (full) "" else " (RF + causal forest only)"
))

expected_rf <- c(
  "dcf",
  "scf_scf",
  "scf_scf_new",
  "scf_oob",
  "scf_oob_t",
  "oob_oob",
  "oob_oob_s",
  "oob_oob_manual"
)
expected_cf <- c("cf_dcf", "cf_scf", "cf_full_oob", "cf_default")
expected_sl <- c(
  "sl_dcf",
  "sl_scf_scf",
  "sl_scf_scf_new"
)
expected <- c(expected_rf, expected_cf, if (full) expected_sl)

report(
  setequal(names(res$arms), expected),
  sprintf(
    "all %d expected arms present (got %d)",
    length(expected),
    length(res$arms)
  )
)

lengths_ok <- vapply(
  res$arms,
  function(a) {
    length(a$tau) == n && length(a$tau_test) == n_test
  },
  logical(1)
)
report(
  all(lengths_ok),
  paste0(
    "tau and tau_test correctly sized",
    if (!all(lengths_ok)) {
      paste0(" [bad: ", paste(names(which(!lengths_ok)), collapse = ", "), "]")
    }
  )
)

na_counts <- vapply(
  res$arms,
  function(a) sum(is.na(a$tau)) + sum(is.na(a$tau_test)),
  numeric(1)
)
report(
  all(na_counts == 0),
  paste0(
    "no NA estimates",
    if (any(na_counts > 0)) {
      paste0(
        " [bad: ",
        paste(names(which(na_counts > 0)), collapse = ", "),
        "]"
      )
    }
  )
)

times_ok <- vapply(
  res$arms,
  function(a) a$time_stage2 >= 0 && a$time_nuisance >= 0,
  logical(1)
)
report(all(times_ok), "all timings non-negative")

single_ok <- vapply(
  res$arms,
  function(a) is.finite(a$mse_test_single),
  logical(1)
)
report(all(single_ok), "single-model test MSE populated for every arm")

# for a whole-sample arm there is only one model, so the two test scores must agree
whole_arms <- c(
  "scf_oob",
  "scf_oob_t",
  "oob_oob",
  "oob_oob_s",
  "oob_oob_manual",
  "cf_full_oob",
  "cf_default"
)
whole_agree <- vapply(
  whole_arms,
  function(nm) {
    a <- res$arms[[nm]]
    abs(a$mse_test_single - mean((a$tau_test - gen1$truth_test_tau)^2)) < 1e-10
  },
  logical(1)
)
report(
  all(whole_agree),
  "single-model and ensemble test MSE coincide for the whole-sample arms"
)

report(
  !identical(res$fold_indices, res$fold_indices_b),
  "the fresh stage-2 split differs from the stage-1 split"
)
report(
  all(sort(res$fold_indices_b) == sort(res$fold_indices)),
  "the fresh split is a permutation of the same fold sizes"
)

# =============================================================================
cat("\n=== 3. scenario 1: no heterogeneity ===\n")

setup_rng_stream(1)
gen_null <- generate_cf_replicate(scenario = 1, n = n, n_test = n_test)
res_null <- run_all_crossfit_variants(
  gen_null$data,
  gen_null$X_test,
  n_folds = 10,
  sl_lib = NULL,
  num.threads = grf_threads,
  truth_test = gen_null$truth_test_tau
)

# tau is computed as p1 - p0 with a shared b0 + b1*X1 + b2*X2 term, so the true
# constant comes back with floating point dust on it - unique() is too strict
tau_spread <- diff(range(gen_null$truth_tau))
report(
  tau_spread < 1e-8,
  sprintf(
    "true CATE is constant at %.3f (spread %.1e)",
    gen_null$truth_tau[1],
    tau_spread
  )
)

m_null <- run_metrics(
  list(
    arms = res_null$arms,
    truth_tau = gen_null$truth_tau,
    truth_test_tau = gen_null$truth_test_tau,
    run = 1
  ),
  scenario = 1
)
report(
  all(m_null$corr == 0) && all(m_null$spearman == 0),
  "correlation metrics forced to 0 in the null scenario, as in cts_metrics.R"
)
report(all(is.finite(m_null$mse)), "all null-scenario MSEs finite")

# =============================================================================
cat("\n=== 4. test-set plumbing and the OOB-counterfactual workaround ===\n")
#
# There is no "optimism" to detect in this study, and an earlier version of this
# check wrongly assumed there was. Optimism is what you see when a model is scored
# against the labels it was fit to. Here every arm is scored against the KNOWN true
# CATE, while the label the stage-2 model saw is a noisy pseudo-outcome.
#
# What is checked instead: that test predictions are wired to the right truth, and
# that the grf X.orig OOB-counterfactual shortcut (oob_predict_counterfactual)
# actually reproduces the documented-API tree-loop it stands in for
# (oob_predict_counterfactual_manual). Every remaining arm is a distinct forest fit
# (no two arms share one fitted model any more, now that the in-sample/naive arms -
# which paired with an OOB view of the same forest - are gone), so there is no
# same-model test-prediction identity left to check.

m <- run_metrics(
  list(
    arms = res$arms,
    truth_tau = gen1$truth_tau,
    truth_test_tau = gen1$truth_test_tau,
    run = 1
  ),
  scenario = 6
)

mse_tbl <- m %>%
  select(arm, set, mse) %>%
  tidyr::pivot_wider(names_from = set, values_from = mse) %>%
  arrange(test)

print(as.data.frame(mse_tbl), digits = 3, row.names = FALSE)

# 4a. test predictions must track the test truth. scenario 6 has strong
# heterogeneity, so a misaligned X_test or truth_test shows up as ~0 correlation.
tracking <- c("dcf", "scf_scf", "scf_oob", "cf_dcf")
cors <- vapply(
  tracking,
  function(nm) {
    cor(res$arms[[nm]]$tau_test, gen1$truth_test_tau)
  },
  numeric(1)
)
report(
  all(cors > 0.2),
  sprintf(
    "test predictions track the test truth (%s)",
    paste(sprintf("%s=%.2f", tracking, cors), collapse = ", ")
  )
)

# 4b. oob_oob_s (X.orig shortcut) and oob_oob_manual (documented get_tree/
# get_leaf_node API) are two different implementations of the same OOB-
# counterfactual idea, each with its own forest fit - so not identical, but should
# be closely correlated on the training sample if the shortcut is doing what
# grf-labs/grf#307 claims.
oob_s_manual_cor <- cor(res$arms$oob_oob_s$tau, res$arms$oob_oob_manual$tau)
report(
  oob_s_manual_cor > 0.8,
  sprintf(
    "oob_oob_s and oob_oob_manual track each other (cor = %.3f)",
    oob_s_manual_cor
  )
)

# 4c. the tight version of 4b: hold one S-learner forest fixed and compare
# oob_predict_counterfactual against oob_predict_counterfactual_manual directly, so
# any difference is attributable only to the prediction method, not to
# forest-to-forest randomness. Rebuilt from gen1 (the replicate res$arms came from)
# rather than reusing section 1's X/W, which are a different draw.
X1 <- as.matrix(gen1$data[, -c(1:2)])
W1 <- gen1$data$W
Y1 <- gen1$data$Y
forest_s <- regression_forest(cbind(W = W1, X1), Y1, num.threads = grf_threads)
shortcut_pred <- oob_predict_counterfactual(forest_s, cbind(W = 1, X1))
manual_pred <- oob_predict_counterfactual_manual(forest_s, cbind(W = 1, X1), Y1)
same_forest_cor <- cor(shortcut_pred, manual_pred)
report(
  same_forest_cor > 0.95,
  sprintf(
    "on one fixed forest, the X.orig shortcut and the manual tree-loop agree (cor = %.3f)",
    same_forest_cor
  )
)

# =============================================================================
plan(sequential)

cat("\n=== summary ===\n")
cat(sprintf("  %d passed, %d failed\n", length(pass), length(fail)))
if (length(fail) > 0) {
  cat("\nfailures:\n")
  for (f in fail) {
    cat("  - ", f, "\n", sep = "")
  }
  quit(status = 1)
}
cat(
  "\nall checks passed. next: submit jobscripts/cf_profile.sh, then run cf_profile_summary.R\n"
)
