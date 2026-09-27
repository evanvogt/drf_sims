##########
# title: verification checks for the crossfitting comparison
##########
# Run this before submitting anything to the HPC:
#
#   Rscript crossfitting/cf_testing.R       # quick: structure + regression check
#   Rscript crossfitting/cf_testing.R full  # adds the SuperLearner family (slow)
#
# Checks, in order:
#   1. the oob_oob arm reproduces R/cate_models.R's dr_random_forest, and
#      propensity trimming is a no-op for the double-crossfit nuisances
#   1b. the per-arm outcome models never see the other arm's outcomes
#   2. every arm returns complete tau / tau_test of the right length
#   3. scenario 1 (no heterogeneity) behaves
#   4. test-set plumbing - test predictions track the test truth. NOT an
#      optimism check: scoring is against the known true CATE, not the labels
#      the models were fit to, so there is no optimism to detect.

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
cat("\n=== 1. regression check: oob_oob against R/cate_models.R ===\n")
# oob_oob is production's dr_random_forest: nuisance_rf (which nuisance_oob_rf
# calls straight through) feeding a whole-sample OOB stage-2 forest. The stage-2
# forest is this folder's own stage2_whole_rf, which masks production's (a
# different signature) once cf_models.R is sourced - so production's is sourced
# into an environment of its own and the two pipelines are compared from the
# same RNG state. They should be bit-identical.

prod <- new.env()
sys.source(here("R", "cate_models.R"), envir = prod)

setup_rng_stream(7)
gen <- generate_cf_replicate(scenario = 8, n = 500, n_test = 1000)
X <- as.matrix(gen$data[, -c(1:2)])
Y <- gen$data$Y
W <- gen$data$W

fold_indices <- sort(seq(nrow(X)) %% 10) + 1
fold_list <- unique(fold_indices)
fold_pairs <- utils::combn(fold_list, 2, simplify = FALSE)

setup_rng_stream(7)
nz_prod <- prod$nuisance_rf(X, Y, W, num.threads = grf_threads)
tau_prod <- prod$stage2_whole_rf(X, nz_prod$po, num.threads = grf_threads)$tau

setup_rng_stream(7)
nz_oob <- nuisance_oob_rf(X, Y, W, num.threads = grf_threads)
tau_oob <- stage2_whole_rf(X, nz_oob$po, gen$X_test, num.threads = grf_threads)$tau_oob

report(identical(nz_prod$po, nz_oob$po),
       "oob_oob pseudo-outcomes are production nuisance_rf's")
report(identical(tau_prod, tau_oob),
       "oob_oob CATE estimates match R/cate_models.R's dr_random_forest")

# trim_ps clamps to [0.05, 0.95]; with W ~ Bernoulli(0.5) no crossfit
# propensity should get anywhere near, so no value may sit on a bound
setup_rng_stream(7)
nz_double <- nuisance_double_rf(X, Y, W, fold_indices, fold_pairs, grf_threads)
ps <- nz_double$W.hat_matrix[!is.na(nz_double$W.hat_matrix)]
cat(sprintf("  double-crossfit propensity range %.3f to %.3f\n", min(ps), max(ps)))
report(all(ps > 0.05 & ps < 0.95),
       "RF propensities lie inside the trimming bounds, so trim_ps changes nothing")

# =============================================================================
cat("\n=== 1b. per-arm outcome models never see the other arm ===\n")
# the crossfit arms' outcome model is t_learner_rf_split: one forest per arm,
# control first. Perturbing the treated units' outcomes must leave the control
# forest - and so Y0.hat - bit-identical, since the control forest is fit from
# the same RNG state on the same rows either way.

in_train <- fold_indices != 1
Y_pert <- Y
Y_pert[W == 1] <- Y_pert[W == 1] + 100

setup_rng_stream(7)
mu <- t_learner_rf_split(X, Y, W, in_train, !in_train, num.threads = grf_threads)
setup_rng_stream(7)
mu_pert <- t_learner_rf_split(X, Y_pert, W, in_train, !in_train, num.threads = grf_threads)

report(identical(mu$Y0.hat, mu_pert$Y0.hat),
       "control-arm predictions unchanged when treated outcomes are perturbed")
report(min(mu_pert$Y1.hat - mu$Y1.hat) > 50,
       "treated-arm predictions move with the treated outcomes")

# =============================================================================
cat("\n=== 2. structure: every arm complete and correctly sized ===\n")

n <- 500
n_test <- 1000
sl_lib <- if (full) {
  sl_libraries(n)
} else {
  NULL
}

setup_rng_stream(1)
gen1 <- generate_cf_replicate(scenario = 8, n = n, n_test = n_test)

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
  "scf_oob",
  "oob_oob"
)
expected_cf <- c("cf_dcf", "cf_scf", "cf_full_oob", "cf_default")
expected_sl <- c(
  "sl_dcf",
  "sl_scf_scf"
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
  "oob_oob",
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
cat("\n=== 4. test-set plumbing ===\n")
#
# No 'optimism' to detect in this study (fixed test-design flaw).
# Test predictions are wired to the right truth.

m <- run_metrics(
  list(
    arms = res$arms,
    truth_tau = gen1$truth_tau,
    truth_test_tau = gen1$truth_test_tau,
    run = 1
  ),
  scenario = 8
)

mse_tbl <- m %>%
  select(arm, set, mse) %>%
  tidyr::pivot_wider(names_from = set, values_from = mse) %>%
  arrange(test)

print(as.data.frame(mse_tbl), digits = 3, row.names = FALSE)

# test predictions must track the test truth. scenario 8 has strong
# heterogeneity, so a misaligned X_test or truth_test shows up as ~0 correlation.
tracking <- c("dcf", "scf_scf", "scf_oob", "oob_oob", "cf_dcf")
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
  "\nall checks passed. next: submit jobscripts/cf_1.sh\n"
)
