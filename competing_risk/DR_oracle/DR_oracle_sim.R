###################
# title: DR oracle 5000 runs - to be a massive array job
# date started: 04/03/2025
# date finished:
# author: Ellie Van Vogt
###################

set.seed(1998)
# libraries
library(dplyr)
library(furrr)
library(future.apply)
library(grf)
library(GenericML)
library(syrup)

# paths 
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# functions
source("live/scripts/competing_risk/functions/collate_predictions.R")
source("live/scripts/competing_risk/DR_oracle/DR_oracle_run.R")

#args and params
args <- commandArgs(trailingOnly = T)
scenario <- as.character(args[1])
n <- as.numeric(args[2])
sim <- as.numeric(args[3])
n_folds <- 10
B <- 200
workers <- 5

#load in the data
datasets <- readRDS(paste0(c("live/data/competing_risk//", scenario, "_", n, ".RDS"), collapse = ""))
datasets <- lapply(datasets, `[[`, 1) # just want the data not the truth

# pick out the data set to be analysed
data <- datasets[[sim]]
#DR oracle run
t0 <- Sys.time()
DR_oracle_res <- DR_oracle_output(data, n_folds, scenario, B, workers)
t1 <- Sys.time()
saveRDS(DR_oracle_res, paste0("live/results/competing_risk/", scenario, "/", n, "/DR_oracle/res_sim_", sim, ".RDS"))
print(paste0(scenario, "_", n, " ran successfully!"))
print(t1-t0)