###############
# script for running all the CATE models in one run - bin outcome
###############

library(dplyr)
library(furrr)
library(grf)
library(GenericML)
library(SuperLearner)
library(here)

# path
path <- here()

# functions
source(here("binary", "bin_dgms.R"))
source(here("binary", "bin_models.R"))
source(here("utils.R"))

# simulation parameters
i <- as.numeric(commandArgs(trailingOnly = T))

workers <- 2

params <- expand.grid(
  scenario = c(1, 3, 8, 9),
  n = c(100, 250, 500, 1000),
  run = c(1:100),
  stringsAsFactors = F
)

# select parameters for current run
param <- params[i,]
print(param)

scenario <- param$scenario
n <- param$n
run <- param$run

n_folds <- ifelse(n == 100, 4, 10)

sl_lib <- c("SL.glm", "SL.glmnet", "SL.earth", "SL.gam", "SL.mean", "SL.randomForest")

# set up simulation seed
setup_rng_stream(run)

# dataset
gen <- generate_binary_scenario_data(scenario, n)

data <- gen$dataset

fmla_info <- get_binary_oracle_info(scenario, gen$bW)

# Run all CATE methods
metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)


results <- run_all_cate_methods(
  data = data, 
  n_folds = n_folds, 
  sl_lib = sl_lib,
  fmla_info = fmla_info
)

results$data <- data
results$truth <- gen$truth

# Save results
output_dir <- file.path(dirname(path), "results", "binary", paste0("scenario_", scenario), n)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(output_dir, paste0("res_sim_", run, ".RDS")))

print("Simulation completed!")