###############
# script for running all the CATE models in one run - cts outcome
###############

library(dplyr)
library(furrr)
library(grf)
library(GenericML)
library(SuperLearner)


path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# functions
source("/rds/general/project/nihr_drf_simulations/live/scripts/drf_sims/continuous/all_methods/cts_dgms.R")
source("/rds/general/project/nihr_drf_simulations/live/scripts/drf_sims/continuous/all_methods/cts_modelling.R")
source("/rds/general/project/nihr_drf_simulations/live/scripts/drf_sims/utils.R")

# arguments to get scenario and simulation number
args <- commandArgs(trailingOnly = T)
scenario <- as.numeric(args[1])
n <- as.numeric(args[2])
sim <- as.numeric(args[3])
n_folds <- ifelse(n == 100, 4, 10)
B <- 200
workers <- 2


# SuperLearner library
sl_lib <- c("SL.glm", "SL.glmnet", "SL.earth", "SL.gam", "SL.mean", "SL.randomForest")


# set up simulation seed
setup_rng_stream(sim)

# dataset
gen <- generate_continuous_scenario_data(scenario, n)

data <- gen$dataset

fmla_info <- get_continuous_oracle_info(scenario, gen$bW)

# Run all CATE methods
results <- run_all_cate_methods(
  data = data, 
  n_folds = n_folds, 
  workers = workers,
  sl_lib = sl_lib,
  fmla_info = fmla_info
)

results$data <- data
results$truth <- gen$truth

# Save results
output_dir <- paste0("live/results/continuous/scenario_", scenario, "/", n, "/all_methods/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, paste0(output_dir, "res_sim_", sim, ".RDS"))

print(paste0("All methods for scenario ", scenario, "_", n, " sim ", sim, " completed successfully!"))

