##########
# title: verification checks for the model evaluation study
##########
# Run this before submitting anything to the HPC:
#
#   Rscript model_evaluation/me_testing.R       # quick: no real nuisance fitting
#   Rscript model_evaluation/me_testing.R full  # adds real XGBoost-CV/H2O AutoML (slow)
#
# Checks, in order:
#   1. me_per_model()'s scoring matches a direct calc_infl_score()/
#      calc_dr_risk() call on the same synthetic fixture - i.e. moving the
#      scoring loop from inside the old calc_metrics() to
#      compute_metrics()'s outer loop (see me_metrics.R) didn't change any
#      number.
#   2. the DGM/design-matrix plumbing - Y, W first; covariates all numeric
#      (R/dgm_scenarios.R's output needs no factor expansion, unlike the old
#      benchtm-based data).
#   3. scenario 1 (no heterogeneity) behaves, same check continuous/ and
#      crossfitting/ both use for this scenario.
#   4. the real pipeline's structure: all 9 candidate models present and
#      correctly sized (always); in "full" mode, both nuisance pipelines
#      present with exactly cv/whole (no infold) and the expected columns,
#      plus one real call to me_per_model() on real output.

library(dplyr)
library(future)
library(future.apply)
library(ranger)
library(glmnet)
library(SuperLearner)
library(xgboost)
library(h2o)
library(caret)
library(here)

source(here("R", "utils.R"))
source(here("model_evaluation", "me_dgms.R"))
source(here("model_evaluation", "me_utils.R"))
source(here("model_evaluation", "me_models.R"))
source(here("model_evaluation", "me_nuisance.R"))
source(here("model_evaluation", "me_config.R"))
source(here("model_evaluation", "me_metrics.R")) # calc_infl_score/calc_dr_risk/me_per_model

full <- "full" %in% commandArgs(trailingOnly = TRUE)

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
cat("\n=== 1. scoring regression check: me_per_model() against a direct call ===\n")
# me_per_model() moved calc_metrics()'s "all 9 models, one vectorized pass"
# scoring into compute_metrics()'s "one model at a time" iteration. This
# checks that move didn't change any number, on a synthetic fixture that
# needs no real model fits or real nuisance estimation.

set.seed(42)
n_fake <- 50
tau_hat_fake <- rnorm(n_fake)
true_tau_fake <- rnorm(n_fake)
Y_fake <- rnorm(n_fake)
W_fake <- rbinom(n_fake, 1, 0.5)

fake_nuisance_df <- function() {
  data.frame(
    tau_T = rnorm(n_fake),
    pi = runif(n_fake, 0.2, 0.8),
    phi = rnorm(n_fake)
  )
}

nuis_cv <- fake_nuisance_df()
nuis_whole <- fake_nuisance_df()

sim_res_fake <- list(
  data = list(Y = Y_fake, W = W_fake),
  nuisances = list(xgb = list(cv = nuis_cv, whole = nuis_whole))
)

scored <- me_per_model(
  list(tau = tau_hat_fake), true_tau_fake, "test_model", sim_res_fake, NULL
)

checks <- list(
  infl_cv_xgb = calc_infl_score(tau_hat_fake, nuis_cv$tau_T, nuis_cv$pi, Y_fake, W_fake),
  dr_cv_xgb = calc_dr_risk(tau_hat_fake, nuis_cv$phi),
  infl_whole_xgb = calc_infl_score(tau_hat_fake, nuis_whole$tau_T, nuis_whole$pi, Y_fake, W_fake),
  dr_whole_xgb = calc_dr_risk(tau_hat_fake, nuis_whole$phi),
  true_pehe = mean((tau_hat_fake - true_tau_fake)^2)
)

for (col in names(checks)) {
  diff <- abs(scored[[col]] - checks[[col]])
  report(diff < 1e-10, sprintf("%s matches a direct call (diff = %.2e)", col, diff))
}

report(ncol(scored) == 5, sprintf(
  "5 score columns from one pipeline/2 fold-types (got %d)", ncol(scored)
))

# =============================================================================
cat("\n=== 2. DGM / design-matrix plumbing ===\n")

setup_rng_stream(1)
gen1 <- generate_me_scenario_data(scenario = 1, n = 100)
data1 <- gen1$dataset

report(
  identical(names(data1)[1:2], c("Y", "W")),
  "Y, W are the first two columns of the generated dataset"
)

design1 <- prepare_design_matrix(data1)
all_numeric <- all(vapply(design1$X, is.numeric, logical(1)))
report(all_numeric, "design matrix covariates are all numeric - no factor expansion needed")
report(nrow(design1$X) == 100, "design matrix has the right number of rows")
report(!is.null(gen1$truth$tau), "generate_me_scenario_data() returns truth$tau")

# =============================================================================
cat("\n=== 3. scenario 1: no heterogeneity ===\n")
# tau is computed as p1 - p0 with a shared b0 + b1*X1 + b2*X2 term, so the
# true constant comes back with floating point dust on it - same check
# continuous/ and crossfitting/ both use.

tau_spread <- diff(range(gen1$truth$tau))
report(
  tau_spread < 1e-8,
  sprintf("true CATE is constant at %.3f (spread %.1e)", gen1$truth$tau[1], tau_spread)
)

# =============================================================================
cat(sprintf(
  "\n=== 4. structural check on the real pipeline%s ===\n",
  if (full) "" else " (candidate models only - rerun with 'full' for nuisances)"
))

n_test <- 250 # smallest grid n, to keep this check as fast as possible
n_folds_test <- 5L # matches me_analysis.R's rule at n = 250

setup_rng_stream(1)
gen4 <- generate_me_scenario_data(scenario = 4, n = n_test)
design4 <- prepare_design_matrix(gen4$dataset)
Y4 <- design4$Y
W4 <- design4$W
X4 <- design4$X
kfolds4 <- split_folds(Y4, k = n_folds_test)

t0 <- Sys.time()
model_list <- run_all_candidate_models(
  Y4, W4, X4, kfolds4$fold_indices, kfolds4$fold_list, kfolds4$fold_pairs
)
elapsed_models <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  candidate-model fitting took %.1fs\n", elapsed_models))

report(
  setequal(names(model_list), CANDIDATE_MODELS),
  sprintf(
    "all %d candidate models present (got %d)",
    length(CANDIDATE_MODELS), length(model_list)
  )
)

lengths_ok <- vapply(model_list, function(m) length(m$tau) == n_test, logical(1))
report(
  all(lengths_ok),
  paste0(
    "all tau vectors correctly sized",
    if (!all(lengths_ok)) {
      paste0(" [bad: ", paste(names(which(!lengths_ok)), collapse = ", "), "]")
    }
  )
)

na_ok <- vapply(model_list, function(m) sum(is.na(m$tau)) == 0, logical(1))
report(
  all(na_ok),
  paste0(
    "no NA tau estimates",
    if (!all(na_ok)) {
      paste0(" [bad: ", paste(names(which(!na_ok)), collapse = ", "), "]")
    }
  )
)

if (full) {
  model_seed <- sample.int(2^31 - 1, 1)
  t0 <- Sys.time()
  nuisances4 <- run_all_nuisance_pipelines(
    X4, Y4, W4, kfolds4$fold_indices, kfolds4$fold_list,
    n_cores = 2, mem = "4G", model_seed = model_seed
  )
  elapsed_nuis <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  nuisance-evaluation pipelines took %.1fs\n", elapsed_nuis))

  report(
    setequal(names(nuisances4), c("xgb", "automl")),
    "both nuisance pipelines (xgb, automl) present"
  )

  fold_types_ok <- vapply(nuisances4, function(p) setequal(names(p), c("cv", "whole")), logical(1))
  report(all(fold_types_ok), "both pipelines expose exactly cv/whole - no infold")

  expected_cols <- c("mu0_T", "mu1_T", "mu0_DR", "mu1_DR", "pi", "mu_DR", "tau_T", "phi", "phi05")
  cols_ok <- vapply(names(nuisances4), function(p) {
    all(vapply(nuisances4[[p]], function(ft_df) setequal(names(ft_df), expected_cols), logical(1)))
  }, logical(1))
  report(all(cols_ok), "every pipeline/fold-type data.frame has the expected 9 columns")

  scored_real <- me_per_model(
    model_list[[1]], gen4$truth$tau, names(model_list)[1],
    list(data = list(Y = Y4, W = W4), nuisances = nuisances4), NULL
  )
  report(
    ncol(scored_real) == 9,
    sprintf("me_per_model() produces 9 score columns on real pipeline output (got %d)", ncol(scored_real))
  )
} else {
  cat("  (nuisance pipelines skipped - rerun with 'full' to exercise real XGBoost-CV/H2O AutoML)\n")
}

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
cat(sprintf(
  "\nall checks passed.%s next: submit jobscripts/me_profile.sh, then run me_profile_summary.R\n",
  if (!full) " run with 'full' to also exercise the real nuisance pipelines, then" else ""
))
