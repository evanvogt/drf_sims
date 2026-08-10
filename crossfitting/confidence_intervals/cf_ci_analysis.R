##########
# title: half-sample bootstrap CI pilot - crossfit-structured RF arms
##########
# Adds confidence intervals to all 12 non-SuperLearner arms of the crossfitting
# comparison. Two structural families, and the bootstrap differs between them:
#
#   5 crossfit-structured arms - dcf, scf_scf, scf_scf_new (family dr_rf) and
#     cf_dcf, cf_scf (family causal_forest). Stage 2 is a genuine per-fold
#     crossfit, so rf_half_boot / cf_half_boot refit per fold against a
#     fold-stratified half sample.
#   7 whole-sample/OOB arms - scf_oob, scf_oob_t, oob_oob, oob_oob_s,
#     oob_oob_manual (dr_rf) and cf_full_oob, cf_default (causal_forest). No
#     fold structure to refit against, so rf_oob_half_boot / cf_oob_half_boot
#     refit one forest per draw on an unstratified half sample. Those arms get
#     two bands out of the same refits (all units vs out-of-half only) plus
#     grf's own OOB variance interval - see crossfitting/README.md.
#
# The 3 SuperLearner arms stay out of scope (not RF-based), which is why this
# calls run_all_crossfit_variants with sl_lib = NULL.
#
# A pilot, not the production run: reduced grid (3 scenarios x 50 runs),
# CI_sf fixed at 0.5 rather than swept. Results land in a brand-new
# ../results/crossfitting_ci/ tree - the production run's
# ../results/crossfitting/ is never read or written here.

library(dplyr)
library(furrr)
library(grf)
library(here)

# path
path <- here()

# functions
source(here("crossfitting", "cf_models.R"))  # generate_cf_replicate, run_all_crossfit_variants
source(here("R", "bootstrap_ci.R"))          # {rf,cf}_half_boot, {rf,cf}_oob_half_boot

# simulation parameters
args <- as.numeric(commandArgs(trailingOnly = T))
i <- args[1]
# trailing args follow cf_analysis.R's convention (sane defaults so a bare
# `Rscript cf_ci_analysis.R 1` still works as a local smoke test); CI_boot is
# an extra knob here, overridable for a fast local run without touching the file
CI_boot <- if (length(args) >= 2 && !is.na(args[2])) args[2] else 200
workers <- if (length(args) >= 3 && !is.na(args[3])) args[3] else 2
grf_threads <- if (length(args) >= 4 && !is.na(args[4])) args[4] else 1
n_test <- 2000

params <- expand.grid(
  scenario = c(1, 6, 9),
  n = 500,
  run = c(1:50),
  stringsAsFactors = F
)

# select parameters for current run
param <- params[i,]
print(param)

scenario <- param$scenario
n <- param$n
run <- param$run

n_folds <- 10L
alpha <- 0.05
CI_sf <- 0.5  # grf's default sample.fraction - fixed for this pilot, no sweep

# set up simulation seed
setup_rng_stream(run)

# data generation - training sample plus an independent test sample from the
# same DGP, drawn in the same stream so the whole replicate is reproducible
gen <- generate_cf_replicate(scenario, n, n_test)

Sys.setenv(OMP_NUM_THREADS = grf_threads)

metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)

# point estimates for all 12 RF/CF arms. sl_lib = NULL drops the SuperLearner
# family, and nothing else about the call differs from cf_analysis.R's - so under
# the same setup_rng_stream(run) seed these arms are bit-identical to the
# production study's, and cf_ci_testing.R check 1 asserts exactly that.
structured <- run_all_crossfit_variants(
  data = gen$data,
  X_test = gen$X_test,
  n_folds = n_folds,
  sl_lib = NULL,
  num.threads = grf_threads,
  truth_test = gen$truth_test_tau  # scoring only, never seen by a model
)

X <- as.matrix(gen$data[, -c(1:2)])
Y <- gen$data$Y
W <- gen$data$W
nz <- structured$nuisances
fold_indices <- structured$fold_indices
fold_list <- unique(fold_indices)
fold_indices_b <- structured$fold_indices_b
fold_list_b <- unique(fold_indices_b)

# table-driven bootstrap wiring: all four bootstraps in R/bootstrap_ci.R share
# the same argument order/count - rf_half_boot(X,Y,W,po,...) ignores Y and W, and
# the two OOB bootstraps ignore the trailing fold arguments - so one loop drives
# every arm.
#
# scf_scf_new is the ONE crossfit row that swaps fold_indices for fold_indices_b.
# It must: its point estimate (stage2_crossfit_rf inside run_all_crossfit_variants)
# was fit under that split, even though its po vector comes from the same
# leave-one-fold-out nuisance (nz_single) as scf_scf.
#
# The OOB rows pass fi = NULL, fl = NULL because a whole-sample arm has no folds.
# cf_default is the one arm with no nuisance stage of its own - nz_cf_default is
# the Y.hat/W.hat grf cross-fit internally, kept by run_all_crossfit_variants so
# this arm's bootstrap holds nuisances fixed like every other one.
boot_spec <- list(
  # per-fold crossfit stage 2
  dcf            = list(fn = rf_half_boot,     arg = nz$nz_double$po,      fi = fold_indices,   fl = fold_list),
  scf_scf        = list(fn = rf_half_boot,     arg = nz$nz_single$po,      fi = fold_indices,   fl = fold_list),
  scf_scf_new    = list(fn = rf_half_boot,     arg = nz$nz_single$po,      fi = fold_indices_b, fl = fold_list_b),
  cf_dcf         = list(fn = cf_half_boot,     arg = nz$nz_double,         fi = fold_indices,   fl = fold_list),
  cf_scf         = list(fn = cf_half_boot,     arg = nz$nz_single,         fi = fold_indices,   fl = fold_list),
  # whole-sample stage 2, OOB predictions
  scf_oob        = list(fn = rf_oob_half_boot, arg = nz$nz_single$po,      fi = NULL,           fl = NULL),
  scf_oob_t      = list(fn = rf_oob_half_boot, arg = nz$nz_single_t$po,    fi = NULL,           fl = NULL),
  oob_oob        = list(fn = rf_oob_half_boot, arg = nz$nz_oob$po,         fi = NULL,           fl = NULL),
  oob_oob_s      = list(fn = rf_oob_half_boot, arg = nz$nz_oob_s$po,       fi = NULL,           fl = NULL),
  oob_oob_manual = list(fn = rf_oob_half_boot, arg = nz$nz_oob_manual$po,  fi = NULL,           fl = NULL),
  cf_full_oob    = list(fn = cf_oob_half_boot, arg = nz$nz_single,         fi = NULL,           fl = NULL),
  cf_default     = list(fn = cf_oob_half_boot, arg = nz$nz_cf_default,     fi = NULL,           fl = NULL)
)

for (nm in names(boot_spec)) {
  spec <- boot_spec[[nm]]
  b <- timed(spec$fn(X, Y, W, spec$arg, structured$arms[[nm]]$tau,
                     CI_boot, CI_sf, alpha, spec$fi, spec$fl))
  # copy every returned band field rather than naming hb_lb/hb_ub, so the OOB
  # arms' extra hb_out_lb/hb_out_ub ride along without a second branch
  for (f in setdiff(names(b$value), "draws")) structured$arms[[nm]][[f]] <- b$value[[f]]
  structured$arms[[nm]]$time_boot <- b$time
  # draws (n x CI_boot) intentionally dropped - not saved, see README note on
  # per-run files staying small; nothing downstream needs them (no MI pooling
  # here, unlike confidence_intervals/missing/ci_example)
}

results <- list(
  arms = structured$arms,
  fold_indices = fold_indices,
  fold_indices_b = fold_indices_b,
  truth_tau = gen$truth_tau,
  truth_test_tau = gen$truth_test_tau,
  bW = gen$bW,
  run = run,
  CI_boot = CI_boot,
  CI_sf = CI_sf,
  alpha = alpha
)

# the data and the nuisance matrices are deliberately not saved, same as
# cf_analysis.R - the replicate is reproducible from run via setup_rng_stream

# Save results
output_dir <- file.path(dirname(path), "results", "crossfitting_ci", paste0("scenario_", scenario))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(output_dir, paste0("res_sim_", run, ".RDS")))

print("Simulation completed!")
