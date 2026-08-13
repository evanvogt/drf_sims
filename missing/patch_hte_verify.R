##########
# title: prove R/patch_hte_tests.R is equivalent to having re-run the study
##########
# Run this before letting the patch touch ../results/missing/. It is the
# evidence behind the claim in missing/binary/README.md that the finished
# results did not need re-running, and it is cheap enough to redo whenever
# R/cate_models.R's test functions change.
#
# Standalone manual tool - never source() this from anything else. run_once()
# below flips the global PROFILES$missing$dr_rf_tests with <<- to run the same
# row both ways; that mutation is reset before the script exits, but only a
# script run start-to-finish gets that guarantee.
#
#   cd missing
#   Rscript patch_hte_verify.R                    # continuous, grid row 1
#   Rscript patch_hte_verify.R binary 3452        # binary, grid row 3452
#   Rscript patch_hte_verify.R continuous 1 nosl  # skip the SuperLearner arm
#
# The argument is a row of study$grid - the same index bin_miss_1.sh /
# cts_miss_1.sh hand to the analysis scripts - so `grid_indices(study, ...)`
# picks out the interesting arms. Worth covering at least:
#   method = "complete_cases"       X is complete, every arm runs
#   method = "none"                 X keeps its NAs, the oracle arms are skipped
#   method = "multiple_imputation"  no nuisances saved, must be skipped cleanly
#
# WHAT IT PROVES
#
# The same simulation is run twice from the same seed, once with
# PROFILES$missing$dr_rf_tests = FALSE (how every existing result file was
# produced) and once TRUE (what the study does now). Then:
#
#  1. every arm's tau is identical between the two runs. This is the claim that
#     matters: the added tests draw no random numbers, so they cannot perturb
#     dr_oracle / dr_semi_oracle / dr_superlearner, which are fitted after
#     dr_random_forest and would otherwise pick up a shifted RNG stream.
#  2. add_hte_tests() applied to the FALSE run reproduces the TRUE run's
#     BLP_whole / independence_cate / independence_po exactly. That is what the
#     patch does to a saved file, so exact agreement means a patched file and a
#     re-run file are the same file - bar dr_random_forest$variance, which the
#     FALSE branch discarded and which nothing reads.
#  3. dr_random_forest$independence_po equals causal_forest$independence_po.
#     Same test, same po, same X - so a mismatch means the patch rebuilt X
#     wrongly from the saved data frame.
#  4. a patched file whose BLP is degenerate (constant tau; run_blp_whole()
#     returns NULL by design - GenericML::BLP() cannot identify beta.2) is
#     still recognised as already_patched. This is a regression check, not a
#     property of this row's own data: it forces tau constant on a COPY of the
#     fitted result and re-derives patch_status_of() from that, so it exercises
#     the failure mode regardless of whether this row's real fit was
#     degenerate. Assigning NULL into a list removes the element, so a naive
#     is.null(BLP_whole) sentinel cannot tell "patched, degenerate" from
#     "never patched" - see the comment above already_patched in
#     R/patch_hte_tests.R for the fix.
#
# `nosl` drops the SuperLearner arm. Use it if the arm crashes locally (see
# "Known issue found while profiling" in missing/binary/README.md); checks 1-3
# still hold across the remaining arms, and dr_oracle / dr_semi_oracle are
# fitted after dr_random_forest too, so check 1 keeps its teeth. Check 4 does
# not touch the SuperLearner arm at all.

library(here)
library(future)

source(here("R", "missingness.R"))      # generate_and_process_data, get_oracle_info
source(here("R", "cate_models.R"))      # cate_methods, PROFILES
source(here("R", "metrics.R"))          # CATE_MODELS
source(here("R", "patch_hte_tests.R"))  # add_hte_tests, patch_status_of
source(here("R", "utils.R"))            # setup_rng_stream

args <- commandArgs(trailingOnly = TRUE)
which_study <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "continuous"
grid_row    <- if (length(args) >= 2) as.integer(args[2]) else 1L
use_sl      <- !("nosl" %in% args)

stopifnot(which_study %in% c("continuous", "binary"))

# Everything below this line differs between the two studies; the estimator code
# path does not, so verifying one verifies the other. Both are offered anyway -
# missing/binary is the study that is already 100% run.
cfg <- list(
  continuous = list(config = "missing/continuous/cts_miss_config.R",
                    set    = "continuous_missing",
                    family = gaussian()),
  binary     = list(config = "missing/binary/bin_miss_config.R",
                    set    = "binary_missing",
                    family = binomial())
)[[which_study]]

env <- new.env()
source(here(cfg$config), local = env)
study <- get("study", envir = env)

param <- study$grid[grid_row, ]
cat("study:", which_study, " grid row:", grid_row, "\n")
print(param)

# Both missing-data studies bake the inverse link into the oracle formula
# string, so the model must not apply it again - see missing/binary/bin_miss_models.R.
ORACLE_LINK <- "identity"
sl_lib <- if (use_sl) {
  c("SL.glm", "SL.glmnet", "SL.earth", "SL.gam", "SL.mean", "SL.ranger")
} else NULL

# Sequential, so the comparison cannot turn on how many workers happened to be
# available. The analysis scripts use multisession; that changes how the
# SuperLearner arm is scheduled, not what either run computes, and both runs
# here are scheduled identically.
plan(sequential)

run_once <- function(dr_rf_tests) {
  # cate_methods() reads PROFILES lexically, and everything here is sourced into
  # the global environment, so this is the switch the whole comparison turns on.
  PROFILES$missing$dr_rf_tests <<- dr_rf_tests

  setup_rng_stream(param$run)

  gen <- generate_and_process_data(
    param$scenario, param$n, set = cfg$set, return_truth = TRUE,
    type = param$type, prop = param$prop, mech = param$mechanism,
    method = param$method, n_imp = 50)

  res <- if (param$method == "multiple_imputation") {
    # Not a simulation of the MI arm - the analysis scripts Rubin-combine 50
    # fits and keep no nuisances. One fit is enough to build the shape the patch
    # has to refuse, which is all this branch is for.
    out <- list()
    out$data <- gen$dataset          # a list of 50, exactly as MI saves it
    out$truth <- gen$truth
    out
  } else {
    fmla_info <- get_oracle_info(param$scenario, gen$bW, set = cfg$set)
    out <- cate_methods(
      gen$dataset, n_folds = 10, sl_lib = sl_lib, fmla_info = fmla_info,
      family = cfg$family, oracle_link = ORACLE_LINK,
      ipw = if (param$method == "IPW") gen$ipw else NULL,
      profile = "missing", num.threads = 1)
    out$data <- gen$dataset
    out$truth <- gen$truth
    out
  }
  res
}

cat("\n--- run 1: dr_rf_tests = FALSE (how the saved results were produced)\n")
res_false <- run_once(FALSE)

cat("\n--- run 2: dr_rf_tests = TRUE (what the study does now)\n")
res_true <- run_once(TRUE)

PROFILES$missing$dr_rf_tests <- TRUE   # leave the session in the real state

# ---- results ----------------------------------------------------------------

failures <- character(0)
check <- function(label, ok) {
  cat(if (isTRUE(ok)) "  PASS  " else "  FAIL  ", label, "\n", sep = "")
  if (!isTRUE(ok)) failures <<- c(failures, label)
}

cat("\n=== the arm the patch refuses ===\n")
if (param$method == "multiple_imputation") {
  check("MI run is reported skipped_mi, not patched or errored",
        patch_status_of(res_false) == "skipped_mi")
  cat("\nMI arms carry no nuisances, so there is nothing to verify beyond the",
      "\nrefusal. See the multiple-imputation note in missing/README.md.\n")
} else {
  check("non-MI run is reported to_patch",
        patch_status_of(res_false) == "to_patch")

  cat("\n=== 1. adding the tests perturbs nothing ===\n")
  arms <- intersect(names(res_true), CATE_MODELS)
  cat("  arms fitted:", paste(arms, collapse = ", "), "\n")
  for (m in arms) {
    check(paste0("tau identical between the two runs: ", m),
          identical(res_false[[m]]$tau, res_true[[m]]$tau))
  }

  cat("\n=== 2. the patch reproduces the re-run exactly ===\n")
  patched <- add_hte_tests(res_false)
  for (fld in c("BLP_whole", "independence_cate", "independence_po")) {
    check(paste0("dr_random_forest$", fld, " matches the re-run"),
          isTRUE(all.equal(patched$dr_random_forest[[fld]],
                           res_true$dr_random_forest[[fld]],
                           tolerance = 0)))
  }
  # Without this, a scenario whose BLP is degenerate (bug L: run_blp_whole
  # returns NULL for a constant tau) would pass check 2 on NULL == NULL and
  # prove nothing. Scenario 1 has no heterogeneity, so verify a scenario that
  # does before believing the result.
  check("BLP_whole is an actual coefficient block, not the bug-L NULL",
        !is.null(patched$dr_random_forest$BLP_whole))
  if (is.null(patched$dr_random_forest$BLP_whole)) {
    cat("        (degenerate tau for this row - re-run against a scenario",
        "with heterogeneity)\n")
  }

  cat("\n=== 3. X was rebuilt as cate_methods() had it ===\n")
  check("dr_random_forest$independence_po == causal_forest$independence_po",
        isTRUE(all.equal(patched$dr_random_forest$independence_po,
                         patched$causal_forest$independence_po, tolerance = 0)))

  cat("\n=== 4. idempotency survives a degenerate BLP ===\n")
  # A copy, not this row's own fit: res_false was left untouched by
  # add_hte_tests() above (R copies on modify), so this only borrows its shape.
  degenerate <- res_false
  degenerate$dr_random_forest$tau <- rep(2.5, length(degenerate$dr_random_forest$tau))
  degenerate_patched <- add_hte_tests(degenerate)

  check("forcing a constant tau genuinely triggers the bug-L NULL (test is valid)",
        is.null(degenerate_patched$dr_random_forest$BLP_whole))
  check("a patched, degenerate-tau file is recognised as already_patched",
        patch_status_of(degenerate_patched) == "already_patched")

  cat("\n=== not recoverable, by design ===\n")
  cat("  dr_random_forest$variance present in the re-run: ",
      !is.null(res_true$dr_random_forest$variance), "\n", sep = "")
  cat("  dr_random_forest$variance present in the patch:  ",
      !is.null(patched$dr_random_forest$variance),
      "  (expected FALSE - see R/patch_hte_tests.R)\n", sep = "")
}

cat("\n")
if (length(failures)) {
  stop(length(failures), " check(s) failed:\n  ", paste(failures, collapse = "\n  "))
}
cat("All checks passed. The patch is equivalent to re-running this row.\n")
