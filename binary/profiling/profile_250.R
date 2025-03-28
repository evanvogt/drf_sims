# test how long running all the 250 scenarios would take

set.seed(1998)
# libraries
library(dplyr)
library(furrr)
library(grf)
library(GenericML)
library(syrup)

# paths 
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# functions
source("live/scripts/binary/functions/collate_predictions.R")
source("live/scripts/binary/CF/CF_run.R")
source("live/scripts/binary/DR_RF/DR_RF_run.R")
source("live/scripts/binary/DR_oracle/DR_oracle_run.R")
source("live/scripts/binary/T_RF/T_RF_run.R")

#args and params
n_folds <- 10
n <- 250
B <- 200
workers <- 10


# do all the big models for one simulated dataset
scens <- paste0("scenario_", seq_along(1:11))
profiling <- syrup(
  lapply(scens, function(scenario) {
    datasets <- readRDS(paste0(c("live/data/binary/", scenario, "_", n, ".RDS"), collapse = ""))
    datasets <- lapply(datasets, `[[`, 1) # just want the data not the truth
    
    data <- datasets[[1]]
    CF_res <- CF_output(data, n_folds, scenario)
    DR_RF_res <- DR_RF_output(data, n_folds, scenario, B, workers)
    DR_oracle_res <- DR_oracle_output(data, n_folds, scenario, B, workers)
    T_RF_output <- T_RF_output(data, n_folds, scenario, B, workers)
    print(paste0(scenario, " ", n, " complete"))
  })
)


saveRDS(profiling, "live/data/binary/profiling/run_250_once.RDS")