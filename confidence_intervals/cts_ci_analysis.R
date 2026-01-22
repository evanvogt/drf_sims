##########
# title: half-sample bootstrap estimation - continuous outcome
##########

# libraries
library(here)
library(dplyr)

# path
path <- here()

# functions
source(here("utils.R"))
source(here("confidence_intervals/cts_ci_dgms.R")) # maybe this can just be the cts folder?
source(here("confidence_intervals/cts_ci_models.R"))

# simulation parameters
i <- as.numeric(commandArgs(trailingOnly = T))

CI_boot <- 200
alpha <- 0.05
n_folds <- 10
workers <- 2

params <- expand.grid(
  scenario = c(1:10),
  n = c(500, 1000),
  CI_sf = seq(0.05, 0.5, 0.05),
  run = c(1:100)
)

# select parameters for current run
param <- params[i,]
print(param)
scenario <- as.numeric(param$scenario)
n <- as.numeric(param$n)
CI_sf <- as.numeric(param$CI_sf)
run <- as.numeric(param$run)

# set up simulation seed
setup_rng_stream(run)

# data generation
gen <- generate_continuous_scenario_data(scenario, n)

data <- gen$dataset

fmla_info <- get_continuous_oracle_info(scenario, gen$bW)


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
  alpha = alpha
)
warnings()

results$data <- data
results$truth <- gen$truth

# save the results
output_dir <- file.path(dirname(path), "results/confidence_intervals/continuous", paste0("scenario_", scenario), n, CI_sf)
dir.create(output_dir, recursive = T, showWarnings = F)
saveRDS(results, file.path(output_dir, paste0("res_sim_", run, ".RDS")))

print("Simulation completed!")