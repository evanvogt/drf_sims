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
source(here("R", "utils.R"))
source(here("binary/bin_config.R"))

# simulation parameters
i <- as.numeric(commandArgs(trailingOnly = T))

workers <- 2

# The parameter grid lives in the study config, so this script and the
# check/collect scripts cannot disagree about what index i means.
param <- study$grid[i, ]
print(param)

scenario <- param$scenario
n <- param$n
run <- param$run

n_folds <- dplyr::case_when(n == 100 ~ 4L, n == 250 ~ 5L, TRUE ~ 10L)

sl_lib <- if (n <= 100) {
  c("SL.glm", "SL.glmnet", "SL.gam", "SL.mean")
} else {
  c("SL.glm", "SL.glmnet", "SL.earth", "SL.gam", "SL.mean", "SL.ranger")
}

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
output_dir <- file.path(
  dirname(path),
  "results",
  "binary",
  paste0("scenario_", scenario),
  n
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(output_dir, paste0("res_sim_", run, ".RDS")))

print("Simulation completed!")
