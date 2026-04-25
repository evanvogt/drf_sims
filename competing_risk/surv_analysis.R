##########
# Title: Competing risks CATEs
##########

library(dplyr)
library(furrr)
library(grf)
library(GenericML)
library(SuperLearner)
library(here)

# Path
path <- here()

# Functions
source(here("competing_risk", "surv_dgm.R"))
source(here("competing_risk", "surv_models.R"))
source(here("utils.R"))

# Simulation parameters
i <- as.numeric(commandArgs(trailingOnly = T))

workers <- 2

horizon <- 28

params <- expand.grid(
  scenario = 1:7,
  censoring = c(TRUE, FALSE),
  n = c(500),
  run = 1:100,
  stringsAsFactors = FALSE
)

# Select parameters for current run
param <- params[i,]
print(param)

scenario <- param$scenario
n <- param$n
censoring <- param$censoring
run <- param$run

n_folds <- ifelse(n < 300, 5, 10)
t0 <- Sys.time()
# Set up simulation seed
setup_rng_stream(run)

# Dataset Generation
metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)

gen <- generate_surv_data(
  scenario = scenario,
  n = n,
  censoring = censoring
)

data <- gen$dataset

# Main analysis
results <- all_cate_surv_models(
  data = data,
  n_folds = n_folds,
  horizon = horizon
)
t1 <- Sys.time()
results$data <- data
results$truth <- gen$truth
print(t1-t0)
# save results
output_dir <- file.path(dirname(path), "results", "competing_risk", paste0("scenario_", scenario), n, paste0("censor_", censoring))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(output_dir, paste0("res_sim_", run, ".RDS")))

print("Simulation completed!")