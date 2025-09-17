###############
# script for running all the CATE models in one run - binary outcome with missing data
# 
###############

library(dplyr)
library(furrr)
library(grf)
library(GenericML)
library(SuperLearner)


path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# functions
source("/rds/general/project/nihr_drf_simulations/live/scripts/drf_sims/missing/binary/binary_dgm_missing.R")
source("/rds/general/project/nihr_drf_simulations/live/scripts/drf_sims/binary/all_methods/CATE_model_functions_bin.R")
source("/rds/general/project/nihr_drf_simulations/live/scripts/drf_sims/utils.R")

# arguments to get scenario and simulation number
args <- commandArgs(trailingOnly = T)
scenario <- as.numeric(args[1])
n <- as.numeric(args[2])
sim <- as.numeric(args[3])
miss_type <- as.character(args[4])
miss_prop <- as.numeric(args[5])
miss_method <- as.character(args[6])
B <- 200
workers <- 10


# SuperLearner library
sl_lib <- c("SL.glm", "SL.glmnet", "SL.earth", "SL.gam", "SL.mean", "SL.randomForest")


# set up simulation seed
setup_rng_stream(sim)

# dataset
gen <- generate_and_process_data(
  scenario, n, return_truth = TRUE, seed = NULL, add_missingness = TRUE, miss_type, miss_prop,
  handle_missing = miss_method
)
n_folds <- ifelse(n < 201, 4, 10)


data <- gen$dataset

fmla_info <- get_binary_oracle_info(scenario, gen$bW)

# Run all CATE methods
results <- run_all_cate_methods(
  data = data, 
  n_folds = n_folds, 
  B = B, 
  workers = workers,
  scenario = scenario,
  n = n,
  sl_lib = sl_lib,
  fmla_info = fmla_info
)

results$data <- data
results$truth <- gen$truth

# Save results
output_dir <- paste0("live/results/missing/binary/", miss_type, "/", miss_prop, "/", miss_method, "/scenario_", scenario, "/", n, "/all_methods/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, paste0(output_dir, "res_sim_", sim, ".RDS"))

print(paste0("All methods for scenario ", scenario, "_", n, " sim ", sim, " completed successfully!"))

