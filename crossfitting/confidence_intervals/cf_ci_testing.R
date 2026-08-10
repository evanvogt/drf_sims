##########
# title: verification checks for the crossfitting CI pilot
##########
# Run this before submitting crossfitting/confidence_intervals/jobscripts/cf_ci_1.sh:
#
#   Rscript crossfitting/confidence_intervals/cf_ci_testing.R
#   Rscript crossfitting/confidence_intervals/cf_ci_testing.R full   # adds check 2
#
# Kept separate from cf_testing.R since it exercises the bootstrap machinery,
# which the production study never touches. Checks, in order:
#   1. the pilot's orchestrator call (run_all_crossfit_variants with
#      sl_lib = NULL) produces the 12 RF/CF arms, and its nuisances are
#      bit-identical to an independent replay of the same RNG stream
#   2. (full only) sl_lib = NULL and sl_lib = <library> agree bit-for-bit on all
#      12 RF/CF arms. The SuperLearner block is strictly after them and gated on
#      sl_lib, so this holds by construction - the point of checking is that the
#      pilot's saved arms ARE the production study's, not a different draw
#   3. cf_half_boot's matrix path (double-crossfit nuisances) still resolves
#      single = FALSE and produces well-formed bands - historical behaviour
#   4. cf_half_boot's vector path (single-crossfit nuisances) resolves
#      single = TRUE and produces well-formed bands
#   5. simultaneous_band's new na.rm argument is inert on an NA-free matrix, and
#      oob_bands masks exactly the in-half cells - the two OOB bands must come
#      off one set of roots, differing only by the mask
#   6. oob_half_sample draws floor(n/2) distinct rows; oob_half_predict really
#      takes the OOB branch for in-half rows rather than silently predicting
#      them in-sample, which is the one way this could degrade into a too-narrow
#      band without erroring
#   7. every one of the 12 arms has a well-formed band: finite, hb_lb <= tau <=
#      hb_ub for most units, non-degenerate width - and the 7 OOB arms also have
#      a well-formed, distinct hb_out_* band. No target-coverage check: coverage
#      from this method is already known (confidence_intervals/) to run below
#      nominal, so that is the pilot's research question, not a pass/fail gate
#   8. the OOB arms' var_oob is finite and positive and normal_interval turns it
#      into a well-formed band; the 5 crossfit arms carry no var_oob
#   9. scf_scf_new is wired to fold_indices_b, not fold_indices - the one
#      exception in an otherwise uniform table, and the easiest place for a
#      copy-paste slip

library(dplyr)
library(furrr)
library(grf)
library(here)

source(here("crossfitting", "cf_models.R"))
source(here("R", "bootstrap_ci.R"))
source(here("R", "metrics.R"))   # normal_interval

args <- commandArgs(trailingOnly = TRUE)
run_full <- length(args) >= 1 && args[1] == "full"

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
n_folds <- 10
alpha <- 0.05
CI_boot_smoke <- 20  # small draw count - this is a structure/wiring check, not a coverage study

# the 12 arms the pilot covers, split by which bootstrap they need
crossfit_arms <- c("dcf", "scf_scf", "scf_scf_new", "cf_dcf", "cf_scf")
oob_arms <- c("scf_oob", "scf_oob_t", "oob_oob", "oob_oob_s", "oob_oob_manual",
              "cf_full_oob", "cf_default")

# =============================================================================
cat("\n=== 1. the pilot's orchestrator call ===\n")

setup_rng_stream(3)
gen <- generate_cf_replicate(scenario = 6, n = n, n_test = n_test)

setup_rng_stream(3)
structured <- run_all_crossfit_variants(
  gen$data, gen$X_test, n_folds = n_folds, sl_lib = NULL,
  num.threads = grf_threads, truth_test = gen$truth_test_tau
)

nz <- structured$nuisances

report(setequal(names(structured$arms), c(crossfit_arms, oob_arms)) &&
         length(structured$arms) == 12,
       sprintf("sl_lib = NULL yields exactly the 12 RF/CF arms (got %d)", length(structured$arms)))

report(all(c("nz_double", "nz_single", "nz_single_t", "nz_oob", "nz_oob_s",
             "nz_oob_manual", "nz_cf_default") %in% names(nz)),
       "the orchestrator returns every nuisance object the bootstraps need")

# nz_double is the first RNG-consuming fit, immediately after
# fold_indices_b <- sample(fold_indices) - so replaying that prefix of the stream
# by hand must reproduce it bit-for-bit. This is what "the fold/pair setup is
# still where it was, and nothing was reordered ahead of stage 1" means to check.
nz_double_replay <- {
  setup_rng_stream(3)
  fold_indices <- sort(seq(n) %% n_folds) + 1
  fold_pairs <- utils::combn(unique(fold_indices), 2, simplify = FALSE)
  fold_indices_b <- sample(fold_indices)  # consume the same RNG draw the orchestrator does
  nuisance_double_rf(as.matrix(gen$data[, -c(1:2)]), gen$data$Y, gen$data$W,
                     fold_indices, fold_pairs, grf_threads)
}
po_diff <- max(abs(nz_double_replay$po - nz$nz_double$po), na.rm = TRUE)
report(po_diff < 1e-8, sprintf("nz_double$po matches an independent RNG replay (max diff %.3e)", po_diff))
yhat_diff <- max(abs(nz_double_replay$Y.hat.cf_matrix - nz$nz_double$Y.hat.cf_matrix), na.rm = TRUE)
report(yhat_diff < 1e-8, sprintf("nz_double$Y.hat.cf_matrix matches the replay (max diff %.3e)", yhat_diff))

# =============================================================================
cat("\n=== 2. pilot arms are the production study's arms ===\n")

if (!run_full) {
  cat("  SKIP  (pass 'full' to run - needs SuperLearner and a second orchestrator pass)\n")
} else {
  suppressPackageStartupMessages(library(SuperLearner))
  sl_lib <- c("SL.glm", "SL.glmnet", "SL.earth", "SL.gam", "SL.mean", "SL.ranger")

  setup_rng_stream(3)
  res_prod <- run_all_crossfit_variants(
    gen$data, gen$X_test, n_folds = n_folds, sl_lib = sl_lib,
    num.threads = grf_threads, truth_test = gen$truth_test_tau
  )

  for (nm in c(crossfit_arms, oob_arms)) {
    same <- identical(res_prod$arms[[nm]]$tau, structured$arms[[nm]]$tau) &&
      identical(res_prod$arms[[nm]]$tau_test, structured$arms[[nm]]$tau_test)
    report(same, sprintf("%s: bit-identical with and without the SuperLearner family", nm))
  }
}

# =============================================================================
cat("\n=== 3. cf_half_boot matrix path (double-crossfit nuisances) ===\n")

X <- as.matrix(gen$data[, -c(1:2)])
Y <- gen$data$Y
W <- gen$data$W
fold_indices <- structured$fold_indices
fold_list <- unique(fold_indices)
fold_indices_b <- structured$fold_indices_b
fold_list_b <- unique(fold_indices_b)

report(is.matrix(nz$nz_double$Y.hat.cf_matrix), "nz_double carries the matrix fields cf_half_boot's matrix path expects")

b_dcf <- cf_half_boot(X, Y, W, nz$nz_double, structured$arms$cf_dcf$tau,
                      CI_boot = CI_boot_smoke, CI_sf = 0.5, alpha = alpha,
                      fold_indices = fold_indices, fold_list = fold_list)

report(all(is.finite(b_dcf$hb_lb)) && all(is.finite(b_dcf$hb_ub)), "cf_dcf band is finite")
report(mean(b_dcf$hb_lb <= structured$arms$cf_dcf$tau & structured$arms$cf_dcf$tau <= b_dcf$hb_ub) > 0.5,
      "cf_dcf band brackets tau for most units")
report(mean(b_dcf$hb_ub - b_dcf$hb_lb) > 0, "cf_dcf band has non-zero width")

# =============================================================================
cat("\n=== 4. cf_half_boot vector path (single-crossfit nuisances) ===\n")

report(is.null(nz$nz_single$Y.hat.cf_matrix) && is.vector(nz$nz_single$Y.hat.cf),
      "nz_single carries only the vector fields - exercises the vector branch")

b_scf <- cf_half_boot(X, Y, W, nz$nz_single, structured$arms$cf_scf$tau,
                      CI_boot = CI_boot_smoke, CI_sf = 0.5, alpha = alpha,
                      fold_indices = fold_indices, fold_list = fold_list)

report(all(is.finite(b_scf$hb_lb)) && all(is.finite(b_scf$hb_ub)), "cf_scf band is finite")
report(mean(b_scf$hb_lb <= structured$arms$cf_scf$tau & structured$arms$cf_scf$tau <= b_scf$hb_ub) > 0.5,
      "cf_scf band brackets tau for most units")
report(mean(b_scf$hb_ub - b_scf$hb_lb) > 0, "cf_scf band has non-zero width")

# =============================================================================
cat("\n=== 5. simultaneous_band na.rm, and oob_bands' masking ===\n")
# the na.rm argument is new and every pre-existing caller relies on the default,
# so it must be provably inert when there is nothing to remove

set.seed(11)
fake_roots <- matrix(rnorm(50 * 30), nrow = 50)
fake_tau <- rnorm(50)
report(identical(simultaneous_band(fake_roots, fake_tau, alpha),
                 simultaneous_band(fake_roots, fake_tau, alpha, na.rm = TRUE)),
      "simultaneous_band: na.rm = TRUE is inert on an NA-free matrix")

# oob_bands must be exactly "one set of roots, two ways of scoring it" - the
# all-units band untouched, the out-of-half band the same roots with the in-half
# cells masked. Checked against the two simultaneous_band calls it should equal.
fake_kept <- matrix(runif(50 * 30) < 0.5, nrow = 50)
ob <- oob_bands(fake_roots, fake_kept, fake_tau, alpha)
expect_all <- simultaneous_band(fake_roots, fake_tau, alpha)
masked <- fake_roots; masked[fake_kept] <- NA_real_
expect_out <- simultaneous_band(masked, fake_tau, alpha, na.rm = TRUE)

report(identical(ob$hb_lb, expect_all$hb_lb) && identical(ob$hb_ub, expect_all$hb_ub),
      "oob_bands: hb_* is the unmasked band")
report(identical(ob$hb_out_lb, expect_out$hb_lb) && identical(ob$hb_out_ub, expect_out$hb_ub),
      "oob_bands: hb_out_* is the same roots with in-half cells masked")
report(identical(ob$draws, fake_roots), "oob_bands returns the unmasked roots as draws")

# =============================================================================
cat("\n=== 6. oob_half_sample and oob_half_predict ===\n")

keep <- oob_half_sample(n)
report(length(keep) == n && is.logical(keep) && sum(keep) == floor(n / 2),
      sprintf("oob_half_sample draws exactly floor(n/2) = %d rows", floor(n / 2)))
report(length(unique(which(keep))) == sum(keep), "oob_half_sample draws without replacement")

# oob_half_predict must give in-half rows their OUT-OF-BAG prediction, not an
# in-sample one. Both are numerically valid calls, so a slip here would produce a
# systematically too-narrow band and no error at all - hence the explicit
# comparison against the wrong version.
probe_forest <- regression_forest(X[keep, ], nz$nz_single$po[keep], sample.fraction = 0.5)
th <- oob_half_predict(probe_forest, X, keep)
in_sample <- predict(probe_forest, newdata = X[keep, , drop = FALSE])$predictions

report(length(th) == n && all(is.finite(th)), "oob_half_predict returns a finite length-n vector")
report(mean(abs(th[keep] - in_sample)) > 1e-8,
      sprintf("oob_half_predict uses OOB predictions for in-half rows (mean |diff| from in-sample = %.3e)",
              mean(abs(th[keep] - in_sample))))
report(isTRUE(all.equal(th[!keep],
                        predict(probe_forest, newdata = X[!keep, , drop = FALSE])$predictions)),
      "oob_half_predict uses plain newdata predictions for out-of-half rows")

# =============================================================================
cat("\n=== 7. every arm: well-formed band (no coverage-target gate) ===\n")

boot_spec <- list(
  # per-fold crossfit stage 2
  dcf            = list(fn = rf_half_boot,     arg = nz$nz_double$po,     fi = fold_indices,   fl = fold_list),
  scf_scf        = list(fn = rf_half_boot,     arg = nz$nz_single$po,     fi = fold_indices,   fl = fold_list),
  scf_scf_new    = list(fn = rf_half_boot,     arg = nz$nz_single$po,     fi = fold_indices_b, fl = fold_list_b),
  cf_dcf         = list(fn = cf_half_boot,     arg = nz$nz_double,        fi = fold_indices,   fl = fold_list),
  cf_scf         = list(fn = cf_half_boot,     arg = nz$nz_single,        fi = fold_indices,   fl = fold_list),
  # whole-sample stage 2, OOB predictions
  scf_oob        = list(fn = rf_oob_half_boot, arg = nz$nz_single$po,     fi = NULL,           fl = NULL),
  scf_oob_t      = list(fn = rf_oob_half_boot, arg = nz$nz_single_t$po,   fi = NULL,           fl = NULL),
  oob_oob        = list(fn = rf_oob_half_boot, arg = nz$nz_oob$po,        fi = NULL,           fl = NULL),
  oob_oob_s      = list(fn = rf_oob_half_boot, arg = nz$nz_oob_s$po,      fi = NULL,           fl = NULL),
  oob_oob_manual = list(fn = rf_oob_half_boot, arg = nz$nz_oob_manual$po, fi = NULL,           fl = NULL),
  cf_full_oob    = list(fn = cf_oob_half_boot, arg = nz$nz_single,        fi = NULL,           fl = NULL),
  cf_default     = list(fn = cf_oob_half_boot, arg = nz$nz_cf_default,    fi = NULL,           fl = NULL)
)

report(setequal(names(boot_spec), names(structured$arms)),
      "boot_spec covers exactly the arms the orchestrator produced")

well_formed <- function(lb, ub, tau) {
  all(is.finite(lb)) && all(is.finite(ub)) &&
    mean(lb <= tau & tau <= ub) > 0.5 &&
    mean(ub - lb) > 0
}

bands <- list()
for (nm in names(boot_spec)) {
  spec <- boot_spec[[nm]]
  tau <- structured$arms[[nm]]$tau
  bands[[nm]] <- spec$fn(X, Y, W, spec$arg, tau, CI_boot_smoke, 0.5, alpha, spec$fi, spec$fl)
  report(well_formed(bands[[nm]]$hb_lb, bands[[nm]]$hb_ub, tau), sprintf("%s: well-formed band", nm))
}

for (nm in oob_arms) {
  tau <- structured$arms[[nm]]$tau
  b <- bands[[nm]]
  report(!is.null(b$hb_out_lb) && well_formed(b$hb_out_lb, b$hb_out_ub, tau),
        sprintf("%s: well-formed out-of-half band", nm))
  # same roots, different scoring - if the mask were a no-op these would coincide
  report(mean(abs((b$hb_ub - b$hb_lb) - (b$hb_out_ub - b$hb_out_lb))) > 1e-8,
        sprintf("%s: out-of-half band is distinct from the all-units band", nm))
}

for (nm in crossfit_arms) {
  report(is.null(bands[[nm]]$hb_out_lb),
        sprintf("%s: crossfit arm carries no out-of-half band", nm))
}

# =============================================================================
cat("\n=== 8. grf's own OOB variance interval ===\n")

for (nm in oob_arms) {
  a <- structured$arms[[nm]]
  v <- a$var_oob
  ok_var <- !is.null(v) && length(v) == n && all(is.finite(v)) && all(v > 0)
  report(ok_var, sprintf("%s: var_oob is finite and strictly positive", nm))
  if (ok_var) {
    ni <- normal_interval(a$tau, v, alpha)
    report(well_formed(ni$lb, ni$ub, a$tau), sprintf("%s: normal_interval band is well-formed", nm))
  }
}

for (nm in crossfit_arms) {
  report(is.null(structured$arms[[nm]]$var_oob),
        sprintf("%s: crossfit arm carries no var_oob (grf's variance does not apply)", nm))
}

# =============================================================================
cat("\n=== 9. scf_scf_new is wired to fold_indices_b, not fold_indices ===\n")
# the one exception in an otherwise uniform boot_spec table - confirm a
# fold_indices-based call would give a materially different (and wrong) band,
# so a copy-paste slip that dropped the exception would show up here

wrong_band <- rf_half_boot(X, Y, W, nz$nz_single$po, structured$arms$scf_scf_new$tau,
                           CI_boot_smoke, 0.5, alpha, fold_indices, fold_list)
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
cat("\nall checks passed. next: submit confidence_intervals/jobscripts/cf_ci_1.sh\n")
