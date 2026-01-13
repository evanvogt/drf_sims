###############
# script for running all the CATE models in one run - TTE sub-distribution  and cause_specific calculations
###############

library(dplyr)
library(furrr)
library(grf)
library(GenericML)
library(randomForestSRC)


path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# functions
source("/rds/general/project/nihr_drf_simulations/live/scripts/drf_sims/competing_risk/no_censoring/surv_dgm_no_cens.R")
source("/rds/general/project/nihr_drf_simulations/live/scripts/drf_sims/competing_risk/no_censoring/csf_models_no_cens.R")
source("/rds/general/project/nihr_drf_simulations/live/scripts/drf_sims/utils.R")

# arguments to get scenario and simulation number
args <- commandArgs(trailingOnly = T)
scenario <- as.numeric(args[1])
n <- as.numeric(args[2])
sim <- as.numeric(args[3])
n_folds <- ifelse(n==250, 5, 10)
workers <- 4

# set up simulation seed
setup_rng_stream(sim)

# dataset
gen <- generate_csh_data(scenario, n)

data <- gen$dataset

# Run all CATE methods
results <- run_all_cate_methods_survival(
  data = data, 
  n_folds = n_folds, 
  workers = workers,
  horizon = 30
)

results$data <- data
results$truth <- gen$truth

# Save results
output_dir <- file.path("live", "results", "competing_risk", paste0("scenario_", scenario), n, "no_censoring")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(output_dir, paste0("res_sim_", sim, ".RDS")))

print(paste0("All methods for scenario ", scenario, "_", n, " sim ", sim, " completed successfully!"))

