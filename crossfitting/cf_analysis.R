##########
# script for running all the crossfitting variants in one run
##########

library(dplyr)
library(furrr)
library(grf)
library(SuperLearner)
library(here)

# path
path <- here()

# functions
source(here("crossfitting", "cf_models.R"))

# simulation parameters
i <- as.numeric(commandArgs(trailingOnly = T))

workers <- 2
grf_threads <- 1
n_test <- 2000

params <- expand.grid(
  scenario = c(1, 4, 6, 9),
  n = 500,
  run = c(1:500),
  stringsAsFactors = F
)

# select parameters for current run
param <- params[i,]
print(param)

scenario <- param$scenario
n <- param$n
run <- param$run

n_folds <- 10L

sl_lib <- c("SL.glm", "SL.glmnet", "SL.earth", "SL.gam", "SL.mean", "SL.ranger")

# set up simulation seed
setup_rng_stream(run)

# data generation - training sample plus an independent test sample from the
# same DGP, drawn in the same stream so the whole replicate is reproducible
gen <- generate_cf_replicate(scenario, n, n_test)

# Run all crossfitting variants.
# multisession workers are new R processes and inherit this, so it keeps
# SL.ranger and the BLAS in step with grf's num.threads rather than each of the
# workers claiming every allocated core.
Sys.setenv(OMP_NUM_THREADS = grf_threads)

metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)

res <- run_all_crossfit_variants(
  data = gen$data,
  X_test = gen$X_test,
  n_folds = n_folds,
  sl_lib = sl_lib,
  num.threads = grf_threads,
  truth_test = gen$truth_test_tau   # scoring only, never seen by a model
)

results <- list(
  arms = res$arms,
  fold_indices = res$fold_indices,
  fold_indices_b = res$fold_indices_b,
  truth_tau = gen$truth_tau,
  truth_test_tau = gen$truth_test_tau,
  bW = gen$bW,
  run = run
)

# the data and the nuisance matrices are deliberately not saved - the replicate is
# reproducible from run via setup_rng_stream, and keeping the per-run files small
# is what lets the collect job run without the 15gb the other studies need

# Save results
output_dir <- file.path(dirname(path), "results", "crossfitting", paste0("scenario_", scenario))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(output_dir, paste0("res_sim_", run, ".RDS")))

print("Simulation completed!")
