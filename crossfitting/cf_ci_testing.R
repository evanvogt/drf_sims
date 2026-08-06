##########
# title: verification checks for the crossfitting CI pilot
##########
# Run this before submitting crossfitting/jobscripts/cf_ci_1.sh:
#
#   Rscript crossfitting/cf_ci_testing.R
#
# Kept separate from cf_testing.R since it exercises different code (the
# bootstrap machinery and run_crossfit_structured_arms, not the full
# run_all_crossfit_variants). Checks, in order:
#   1. run_crossfit_structured_arms's nz_double/nz_single nuisances are
#      bit-identical to run_all_crossfit_variants's under the same seed
#      (nothing precedes them in the RNG stream in either orchestrator); its
#      dcf/cf_dcf/etc stage-2 estimates are NOT bit-identical (skipping the 4
#      out-of-scope nuisance fits between nz_single and stage 2 shifts every
#      later forest fit's RNG position - see the docstring in cf_models.R),
#      but should still agree closely as estimators of the same target
#   2. cf_half_boot's matrix path (double-crossfit nuisances) still resolves
#      single = FALSE and produces well-formed bands - the historical
#      behaviour must be unchanged
#   3. cf_half_boot's new vector path (single-crossfit nuisances) resolves
#      single = TRUE and produces well-formed bands
#   4. every in-scope arm's band is well-formed: finite, hb_lb <= tau <=
#      hb_ub for most units, non-degenerate width. No target-coverage check -
#      coverage from this method is already known (confidence_intervals/) to
#      run below nominal, so that is the pilot's research question, not a
#      pass/fail gate here.
#   5. scf_scf_new is wired to fold_indices_b, not fold_indices - the one
#      exception in an otherwise uniform table, and the easiest place for a
#      copy-paste slip

library(dplyr)
library(furrr)
library(grf)
library(here)

source(here("crossfitting", "cf_models.R"))
source(here("R", "bootstrap_ci.R"))

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

n <- 500
n_test <- 1000
CI_boot_smoke <- 20  # small draw count - this is a structure/wiring check, not a coverage study

# =============================================================================
cat("\n=== 1. regression check: run_crossfit_structured_arms against run_all_crossfit_variants ===\n")

setup_rng_stream(3)
gen <- generate_cf_replicate(scenario = 6, n = n, n_test = n_test)

setup_rng_stream(3)
res_full <- run_all_crossfit_variants(
  gen$data, gen$X_test, n_folds = 10, sl_lib = NULL,
  num.threads = grf_threads, truth_test = gen$truth_test_tau
)

setup_rng_stream(3)
structured <- run_crossfit_structured_arms(
  gen$data, gen$X_test, n_folds = 10,
  num.threads = grf_threads, truth_test = gen$truth_test_tau
)

report(
  identical(structured$fold_indices, res_full$fold_indices) &&
    identical(structured$fold_indices_b, res_full$fold_indices_b),
  "fold_indices and fold_indices_b match between the two orchestrators"
)

# nz_double/nz_single are the first RNG-consuming fits in BOTH orchestrators,
# immediately after fold_indices_b <- sample(fold_indices) - nothing precedes
# them in either, so they must be bit-identical. This is what "the extraction
# preserved the fold/pair setup and didn't accidentally reorder anything before
# stage 1" actually means to check. Stage-2 outputs (tau, tau_test_folds, ...)
# diverge downstream (see file header), so nuisances are compared directly by
# recomputing them via the same call run_all_crossfit_variants itself makes.
nz_double_full <- {
  setup_rng_stream(3)
  fold_indices <- sort(seq(n) %% 10) + 1
  fold_pairs <- utils::combn(unique(fold_indices), 2, simplify = FALSE)
  fold_indices_b <- sample(fold_indices)  # consume the same RNG draw run_all_crossfit_variants does
  nuisance_double_rf(as.matrix(gen$data[, -c(1:2)]), gen$data$Y, gen$data$W,
                     fold_indices, fold_pairs, grf_threads)
}
po_diff <- max(abs(nz_double_full$po - structured$nz_double$po), na.rm = TRUE)
report(po_diff < 1e-8, sprintf("nz_double$po is bit-identical between the two orchestrators (max diff %.3e)", po_diff))
yhat_diff <- max(abs(nz_double_full$Y.hat.cf_matrix - structured$nz_double$Y.hat.cf_matrix), na.rm = TRUE)
report(yhat_diff < 1e-8, sprintf("nz_double$Y.hat.cf_matrix is bit-identical (max diff %.3e)", yhat_diff))

# dcf/cf_dcf/etc stage-2 estimates diverge from here on (see file header) -
# not a bug, but they should still be highly correlated as estimators of the
# same target on the same data
for (nm in c("dcf", "scf_scf", "scf_scf_new", "cf_dcf", "cf_scf")) {
  agreement <- cor(res_full$arms[[nm]]$tau, structured$arms[[nm]]$tau)
  report(agreement > 0.8,
        sprintf("%s: run_crossfit_structured_arms agrees with run_all_crossfit_variants (cor = %.3f, not bit-identical by design)", nm, agreement))
}

# =============================================================================
cat("\n=== 2. cf_half_boot matrix path (double-crossfit nuisances) ===\n")

X <- as.matrix(gen$data[, -c(1:2)])
Y <- gen$data$Y
W <- gen$data$W
fold_indices <- structured$fold_indices
fold_list <- unique(fold_indices)
fold_indices_b <- structured$fold_indices_b
fold_list_b <- unique(fold_indices_b)

report(is.matrix(structured$nz_double$Y.hat.cf_matrix), "nz_double carries the matrix fields cf_half_boot's matrix path expects")

b_dcf <- cf_half_boot(X, Y, W, structured$nz_double, structured$arms$cf_dcf$tau,
                      CI_boot = CI_boot_smoke, CI_sf = 0.5, alpha = 0.05,
                      fold_indices = fold_indices, fold_list = fold_list)

report(all(is.finite(b_dcf$hb_lb)) && all(is.finite(b_dcf$hb_ub)), "cf_dcf band is finite")
report(mean(b_dcf$hb_lb <= structured$arms$cf_dcf$tau & structured$arms$cf_dcf$tau <= b_dcf$hb_ub) > 0.5,
      "cf_dcf band brackets tau for most units")
report(mean(b_dcf$hb_ub - b_dcf$hb_lb) > 0, "cf_dcf band has non-zero width")

# =============================================================================
cat("\n=== 3. cf_half_boot vector path (single-crossfit nuisances) ===\n")

report(is.null(structured$nz_single$Y.hat.cf_matrix) && is.vector(structured$nz_single$Y.hat.cf),
      "nz_single carries only the vector fields - exercises the new branch")

b_scf <- cf_half_boot(X, Y, W, structured$nz_single, structured$arms$cf_scf$tau,
                      CI_boot = CI_boot_smoke, CI_sf = 0.5, alpha = 0.05,
                      fold_indices = fold_indices, fold_list = fold_list)

report(all(is.finite(b_scf$hb_lb)) && all(is.finite(b_scf$hb_ub)), "cf_scf band is finite")
report(mean(b_scf$hb_lb <= structured$arms$cf_scf$tau & structured$arms$cf_scf$tau <= b_scf$hb_ub) > 0.5,
      "cf_scf band brackets tau for most units")
report(mean(b_scf$hb_ub - b_scf$hb_lb) > 0, "cf_scf band has non-zero width")

# =============================================================================
cat("\n=== 4. every in-scope arm: well-formed band (no coverage-target gate) ===\n")

boot_spec <- list(
  dcf         = list(fn = rf_half_boot, arg = structured$nz_double$po, fi = fold_indices,   fl = fold_list),
  scf_scf     = list(fn = rf_half_boot, arg = structured$nz_single$po, fi = fold_indices,   fl = fold_list),
  scf_scf_new = list(fn = rf_half_boot, arg = structured$nz_single$po, fi = fold_indices_b, fl = fold_list_b),
  cf_dcf      = list(fn = cf_half_boot, arg = structured$nz_double,    fi = fold_indices,   fl = fold_list),
  cf_scf      = list(fn = cf_half_boot, arg = structured$nz_single,    fi = fold_indices,   fl = fold_list)
)

bands <- list()
for (nm in names(boot_spec)) {
  spec <- boot_spec[[nm]]
  bands[[nm]] <- spec$fn(X, Y, W, spec$arg, structured$arms[[nm]]$tau,
                        CI_boot_smoke, 0.5, 0.05, spec$fi, spec$fl)
  ok_finite <- all(is.finite(bands[[nm]]$hb_lb)) && all(is.finite(bands[[nm]]$hb_ub))
  ok_bracket <- mean(bands[[nm]]$hb_lb <= structured$arms[[nm]]$tau &
                       structured$arms[[nm]]$tau <= bands[[nm]]$hb_ub) > 0.5
  ok_width <- mean(bands[[nm]]$hb_ub - bands[[nm]]$hb_lb) > 0
  report(ok_finite && ok_bracket && ok_width, sprintf("%s: well-formed band", nm))
}

# =============================================================================
cat("\n=== 5. scf_scf_new is wired to fold_indices_b, not fold_indices ===\n")
# the one exception in an otherwise uniform boot_spec table - confirm a
# fold_indices-based call would give a materially different (and wrong) band,
# so a copy-paste slip that dropped the exception would show up here

wrong_band <- rf_half_boot(X, Y, W, structured$nz_single$po, structured$arms$scf_scf_new$tau,
                           CI_boot_smoke, 0.5, 0.05, fold_indices, fold_list)
right_band <- bands$scf_scf_new

# fold_indices and fold_indices_b are different permutations of the same fold
# sizes (checked in cf_testing.R), so a half-sample bootstrap stratified by the
# wrong split draws a different set of folds and should give a numerically
# different band width from the correct one
width_diff <- mean(abs((right_band$hb_ub - right_band$hb_lb) - (wrong_band$hb_ub - wrong_band$hb_lb)))
report(width_diff > 1e-6,
      sprintf("scf_scf_new's fold_indices_b band differs from a fold_indices-based band (mean |width diff| = %.3e)", width_diff))

# =============================================================================
plan(sequential)

cat("\n=== summary ===\n")
cat(sprintf("  %d passed, %d failed\n", length(pass), length(fail)))
if (length(fail) > 0) {
  cat("\nfailures:\n")
  for (f in fail) cat("  - ", f, "\n", sep = "")
  quit(status = 1)
}
cat("\nall checks passed. next: submit jobscripts/cf_ci_1.sh\n")
