##########
# title: half-sample bootstrap estimation - continuous outcome
##########

# libraries
library(here)
library(dplyr)

# path
path <- here()

# functions
source(here("R", "utils.R"))
source(here("confidence_intervals", "continuous", "cts_ci_dgms.R")) # maybe this can just be the cts folder?
source(here("confidence_intervals", "continuous", "cts_ci_models.R"))
source(here("confidence_intervals/continuous/cts_ci_config.R"))

# simulation parameters
i <- as.numeric(commandArgs(trailingOnly = T))

CI_boot <- 200
alpha <- 0.05
n_folds <- 10
workers <- 2

# The parameter grid lives in the study config, so this script and the
# check/collect scripts cannot disagree about what index i means.
param <- study$grid[i, ]
print(param)

scenario <- param$scenario
n <- param$n
CI_sf <- param$CI_sf
run <- param$run

# set up simulation seed
setup_rng_stream(run)

# data generation
gen <- generate_continuous_scenario_data(scenario, n)

data <- gen$dataset

fmla_info <- get_continuous_oracle_info(scenario, gen$bW)

# fixed covariate-grid query points, for the grid-based (as opposed to
# per-unit) simultaneous coverage - see R/dgm_scenarios.R::build_query_grid.
# Depends only on scenario (which of X3/X4/X5 are active) and gen$bW, not on
# n/CI_sf/run, so rebuilding it every replicate is cheap and keeps this script
# self-contained, same as fmla_info above.
Z_query <- build_query_grid(scenario, set = "continuous",
                            covariate_names = names(data)[-c(1, 2)])
grid_truth <- build_query_grid_truth(scenario, set = "continuous", gen$bW, Z_query)


# Set up parallelisation
metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)


# run the methods and get CIs
results <- run_all_cate_methods(
  data = data,
  n_folds = n_folds,
  fmla_info = fmla_info,
  CI_boot = CI_boot,
  CI_sf = CI_sf,
  alpha = alpha,
  Z_query = Z_query
)
warnings()

results$data <- data
results$truth <- gen$truth
results$Z_query <- Z_query
results$grid_truth <- grid_truth

# save the results
output_dir <- file.path(dirname(path), "results/confidence_intervals/continuous", paste0("scenario_", scenario), n, CI_sf)
dir.create(output_dir, recursive = T, showWarnings = F)
saveRDS(results, file.path(output_dir, paste0("res_sim_", run, ".RDS")))

print("Simulation completed!")