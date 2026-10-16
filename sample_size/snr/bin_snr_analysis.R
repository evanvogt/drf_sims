###############
# correlated binary SNR check - DR-learners vs causal forest
###############
# ../bin_corr_analysis.R on the binary_snr set (rho = 0.5), with the effect
# modifiers scaled by hte_mult, bW calibrated at SNR_CALIB_N and only the arms
# in SNR_MODELS (bin_snr_config.R).
#   Rscript bin_snr_analysis.R <i> [grf_threads]
# There is no workers argument: only the SuperLearner arm, left out here, runs
# through future_map(), so a future plan would start workers nothing uses.

library(dplyr)
library(grf)
library(here)

# functions
source(here("R", "dgm_scenarios.R"))
source(here("R", "cate_models.R"))
source(here("R", "utils.R"))
source(here("sample_size", "correlated", "binary", "snr", "bin_snr_config.R"))

# simulation parameters
args <- as.numeric(commandArgs(trailingOnly = T))
i <- args[1]
# grf_threads defaults to 1, as bin_corr_analysis.R; the jobscripts pass ncpus
grf_threads <- if (length(args) >= 2 && !is.na(args[2])) args[2] else 1

param <- study$grid[i, ]
print(param)

scenario <- param$scenario
n <- param$n
hte_mult <- param$hte_mult
run <- param$run

# fold_indices only: no arm here crossfits (the SuperLearner arm is left out),
# kept as bin_corr_analysis.R sets it so the saved fold_indices match
n_folds <- dplyr::case_when(n == 100 ~ 4L, n == 250 ~ 5L, TRUE ~ 10L)

t0 <- Sys.time()
# the run alone, as the main study: the n = 500 cells of a run are paired
setup_rng_stream(run)

# dataset - bW calibrated at SNR_CALIB_N, so every n shares one DGM
gen <- generate_scenario_data(scenario, n, set = SNR_SET, calib_n = SNR_CALIB_N,
                              hte_mult = hte_mult)

data <- gen$dataset

fmla_info <- get_oracle_info(scenario, gen$bW, SNR_SET, hte_mult = hte_mult)

# keeps the BLAS in step with grf's num.threads
Sys.setenv(OMP_NUM_THREADS = grf_threads)

# SuperLearner is the last arm cate_methods() runs, so leaving it out does not
# move the RNG stream of the arms before it: at hte_mult = 1 and n = 500 this
# reproduces bin_corr_analysis.R's arms at rho = 0.5 exactly
results <- cate_methods(
  data = data,
  n_folds = n_folds,
  sl_lib = NULL,
  fmla_info = fmla_info,
  family = binomial(),
  profile = "base",
  num.threads = grf_threads,
  models = SNR_MODELS
)
t1 <- Sys.time()

results$data <- data
results$truth <- gen$truth
print(t1 - t0)

# Save results, under the path check/collect expect (combo_dir())
output_dir <- combo_dir(study, param)
save_run(results, output_dir, run)

print("Simulation completed!")
