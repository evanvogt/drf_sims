###################
# title: running cate models across all simulations for one sample size, one sim at a time
# date started: 04/03/2025
# date finished:
# author: Ellie Van Vogt
###################

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
source("live/scripts/functions/collate_predictions.R")
source("live/scripts/CF/CF_run.R")
source("live/scripts/DR_RF/DR_RF_run.R")
source("live/scripts/DR_oracle/DR_oracle_run.R")
source("live/scripts/T_RF/T_RF_run.R")

#args and params
args <- commandArgs(trailingOnly = T)
n <- as.numeric(args[1])
sim <- as.numeric(args[2])

n_folds <- 10
B <- 200
workers <- 10


# do all the big models for one simulated dataset at a time - fix the sample size and go across the scenarios
scens <- paste0("scenario_", seq_along(1:11))

lapply(scens, function(scenario) {
  datasets <- readRDS(paste0(c("live/data/", scenario, "_", n, ".RDS"), collapse = ""))
  datasets <- lapply(datasets, `[[`, 1) # just want the data not the truth
  
  data <- datasets[[sim]]
  CF_res <- CF_output(data, n_folds, scenario)
  saveRDS(CF_res, paste0("live/results/", scenario, "/", n, "/CF/res_sim_", sim, ".RDS"))
  
  DR_RF_res <- DR_RF_output(data, n_folds, scenario, B, workers)
  saveRDS(DR_RF_res, paste0("live/results/", scenario, "/", n, "/DR_RF/res_sim_", sim, ".RDS"))
  
  DR_oracle_res <- DR_oracle_output(data, n_folds, scenario, B, workers)
  saveRDS(DR_oracle_res, paste0("live/results/", scenario, "/", n, "/DR_oracle/res_sim_", sim, ".RDS"))
  
  T_RF_res <- T_RF_output(data, n_folds, scenario, B, workers)
  saveRDS(T_RF_res, paste0("live/results/", scenario, "/", n, "/T_RF/res_sim_", sim, ".RDS"))
})



