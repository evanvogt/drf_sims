##########
# title: half-sample bootstrap CI pilot - crossfit-structured RF arms
##########
# Adds half-sample bootstrap confidence intervals to the 5 arms of the
# crossfitting comparison whose stage 2 is a genuine per-fold crossfit:
# dcf, scf_scf, scf_scf_new (family dr_rf) and cf_dcf, cf_scf (family
# causal_forest). The other 10 arms (6 whole-sample/OOB, no per-fold
# structure to bootstrap the way cf_half_boot/rf_half_boot expect; 3
# SuperLearner, not RF-based) are out of scope - see crossfitting/README.md.
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
source(here("crossfitting", "cf_models.R"))  # generate_cf_replicate, run_crossfit_structured_arms
source(here("R", "bootstrap_ci.R"))          # rf_half_boot, cf_half_boot

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

# point estimates for exactly the 5 crossfit-structured arms
structured <- run_crossfit_structured_arms(
  data = gen$data,
  X_test = gen$X_test,
  n_folds = n_folds,
  num.threads = grf_threads,
  truth_test = gen$truth_test_tau  # scoring only, never seen by a model
)

X <- as.matrix(gen$data[, -c(1:2)])
Y <- gen$data$Y
W <- gen$data$W
fold_indices <- structured$fold_indices
fold_list <- unique(fold_indices)
fold_indices_b <- structured$fold_indices_b
fold_list_b <- unique(fold_indices_b)

# table-driven bootstrap wiring: rf_half_boot(X,Y,W,po,...) and
# cf_half_boot(X,Y,W,nuisances,...) share the same argument order/count, so one
# loop drives both. scf_scf_new is the ONE row that swaps fold_indices for
# fold_indices_b - it must, since its point estimate (stage2_crossfit_rf inside
# run_crossfit_structured_arms) was fit under that split, even though its po
# vector comes from the same leave-one-fold-out nuisance (nz_single) as scf_scf.
boot_spec <- list(
  dcf         = list(fn = rf_half_boot, arg = structured$nz_double$po, fi = fold_indices,   fl = fold_list),
  scf_scf     = list(fn = rf_half_boot, arg = structured$nz_single$po, fi = fold_indices,   fl = fold_list),
  scf_scf_new = list(fn = rf_half_boot, arg = structured$nz_single$po, fi = fold_indices_b, fl = fold_list_b),
  cf_dcf      = list(fn = cf_half_boot, arg = structured$nz_double,    fi = fold_indices,   fl = fold_list),
  cf_scf      = list(fn = cf_half_boot, arg = structured$nz_single,    fi = fold_indices,   fl = fold_list)
)

for (nm in names(boot_spec)) {
  spec <- boot_spec[[nm]]
  b <- timed(spec$fn(X, Y, W, spec$arg, structured$arms[[nm]]$tau,
                     CI_boot, CI_sf, alpha, spec$fi, spec$fl))
  structured$arms[[nm]]$hb_lb <- b$value$hb_lb
  structured$arms[[nm]]$hb_ub <- b$value$hb_ub
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
