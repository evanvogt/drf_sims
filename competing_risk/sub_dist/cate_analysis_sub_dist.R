###############
# script for running all the CATE models in one run - TTE subdistribution effect
# 
###############

library(dplyr)
library(furrr)
library(grf)
library(GenericML)
library(randomForestSRC)


path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# functions
source("/rds/general/project/nihr_drf_simulations/live/scripts/drf_sims/competing_risk/surv_dgm_functions.R")
source("/rds/general/project/nihr_drf_simulations/live/scripts/drf_sims/competing_risk/sub_dist/model_functions_sub_dist.R")
source("/rds/general/project/nihr_drf_simulations/live/scripts/drf_sims/utils.R")

# arguments to get scenario and simulation number
args <- commandArgs(trailingOnly = T)
scenario <- as.numeric(args[1])
n <- as.numeric(args[2])
sim <- as.numeric(args[3])
n_folds <- ifelse(n == 100, 4, 10)
workers <- 10

# set up simulation seed
setup_rng_stream(sim)

# dataset
gen <- generate_survival_scenario_data(scenario, n)

data <- gen$dataset

fmla_info <- get_survival_oracle_info(scenario)

# Run all CATE methods
results <- run_all_cate_methods_survival(
  data = data, 
  n_folds = n_folds, 
  workers = workers,
  horizon = 30,
  fmla_info = fmla_info
)

# Save results
output_dir <- paste0("live/results/competing_risk/sub_dist/scenario_", scenario, "/", n, "/all_methods/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, paste0(output_dir, "res_sim_", sim, ".RDS"))

print(paste0("All methods for scenario ", scenario, "_", n, " sim ", sim, " completed successfully!"))

