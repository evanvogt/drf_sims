###############
# script for running all the CATE models in one run - correlated bin outcome
###############
# The retired binary/bin_analysis.R (tag independent-ss-final) on the
# correlated sets: the same estimators (bin_corr_models.R, formerly
# binary/bin_models.R, family = binomial()), with rho read from the grid row.

library(dplyr)
library(furrr)
library(grf)
library(GenericML)
library(SuperLearner)
library(here)

# functions
source(here("sample_size", "correlated", "binary", "bin_corr_dgms.R"))
source(here("sample_size", "correlated", "binary", "bin_corr_models.R"))
source(here("R", "utils.R"))
source(here("sample_size", "correlated", "binary", "bin_corr_config.R"))

# simulation parameters
args <- as.numeric(commandArgs(trailingOnly = T))
i <- args[1]
# workers/grf_threads default to 2/1 so a bare `Rscript bin_corr_analysis.R 1`
# still works as a local smoke test; the jobscripts supply both explicitly
workers <- if (length(args) >= 2 && !is.na(args[2])) args[2] else 2
grf_threads <- if (length(args) >= 3 && !is.na(args[3])) args[3] else 1

# The parameter grid lives in the study config, so this script and the
# check/collect scripts cannot disagree about what index i means.
param <- study$grid[i, ]
print(param)

scenario <- param$scenario
n <- param$n
run <- param$run
rho <- param$rho

n_folds <- dplyr::case_when(n == 100 ~ 4L, n == 250 ~ 5L, TRUE ~ 10L)

# per-nuisance SuperLearner libraries, smaller at n <= 100 - see R/sl_library.R
sl_lib <- sl_libraries(n)

# set up simulation seed - the run alone, so each rho sees the same draws
setup_rng_stream(run)

# dataset
gen <- generate_binary_scenario_data(scenario, n, rho)

data <- gen$dataset

fmla_info <- get_binary_oracle_info(scenario, gen$bW, rho)

# Run all CATE methods
# multisession workers are new R processes and inherit this, so it keeps
# SL.ranger and the BLAS in step with grf's num.threads - as binary/.
Sys.setenv(OMP_NUM_THREADS = grf_threads)

metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)


results <- run_all_cate_methods(
  data = data,
  n_folds = n_folds,
  sl_lib = sl_lib,
  fmla_info = fmla_info,
  num.threads = grf_threads
)

results$data <- data
results$truth <- gen$truth

# Save results, under the path check/collect expect (combo_dir())
output_dir <- combo_dir(study, param)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(output_dir, paste0("res_sim_", run, ".RDS")))

print("Simulation completed!")
