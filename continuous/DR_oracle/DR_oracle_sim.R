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
library(grf)
library(GenericML)
library(syrup)

# paths 
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# functions
source("live/scripts/functions/collate_predictions.R")
source("live/scripts/DR_oracle/DR_oracle_run.R")

#args and params
args <- commandArgs(trailingOnly = T)
scenario <- as.character(args[1])
n <- as.numeric(args[2])
sim <- as.numeric(args[3])
n_folds <- 10
B <- 200
workers <- 5

#load in the data
datasets <- readRDS(paste0(c("live/data/", scenario, "_", n, ".RDS"), collapse = ""))
datasets <- lapply(datasets, `[[`, 1) # just want the data not the truth

# pick out the data set to be analysed
data <- datasets[[sim]]
#DR oracle run
DR_oracle_res <- DR_oracle_output(data, n_folds, scenario, B, workers)
saveRDS(DR_oracle_res, paste0("live/results/", scenario, "/", n, "/DR_oracle/res_sim_", sim, ".RDS"))
print(paste0(scenario, "_", n, " ran successfully!"))