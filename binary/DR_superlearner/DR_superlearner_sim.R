###################
# title: DR superlearner
# date started: 04/03/2025
# date finished:
# author: Ellie Van Vogt
###################

set.seed(1998)
# libraries
library(dplyr)
library(furrr)
library(SuperLearner)
library(GenericML)

# paths 
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# functions
source("live/scripts/binary/functions/collate_predictions.R")
source("live/scripts/binary/DR_superlearner/DR_superlearner_run.R")

#args and params
args <- commandArgs(trailingOnly = T)
scenario <- as.character(args[1])
n <- as.numeric(args[2])
sim <- as.numeric(args[3])
n_folds <- ifelse(n == "100", 4, 10)
B <- 200
workers <- 5

# load the data
datasets <- readRDS(paste0(c("live/data/binary/", scenario, "_", n, ".RDS"), collapse = ""))
datasets <- lapply(datasets, `[[`, 1) # just want the data not the truth

# pick out the data set to be analysed
data <- datasets[[sim]]

# functions for the superlearner
sl_lib <- c("SL.glm", "SL.glmnet", "SL.earth", "SL.gam", "SL.mean", "SL.randomForest")

#DR run
DR_SL_res <- DR_SL_output(data, n_folds, scenario, B, workers, sl_lib)
saveRDS(DR_SL_res, paste0("live/results/binary/", scenario, "/", n, "/DR_superlearner/res_sim_", sim, ".RDS"))
print(paste0(scenario, "_", n, " ran successfully!"))