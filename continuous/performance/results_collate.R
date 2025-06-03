###############
# title: collecting all the results into one file
# date started: 20/02/25
# date finished:
# author: Ellie Van Vogt
###############

rm(list = ls(all = TRUE))

#libraries
library(dplyr)
library(tidyr)

# paths
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path) 

# functions
source("live/scripts/continuous/functions/extract_results.R")

# headings
# i thjikn results up to 1000 will be available for all the models fingers crossed
samplesizes <- c(250, 500, 1000)#, 5000)
models <- c("CF", "DR_RF", "DR_oracle", "T_RF")#, "Logistic", "H_lasso")
scens <- paste0("scenario_", seq_along(1:10))

# combine results by type I guess so that it's easier to do the performance metrics
cate_preds <- extract_results("tau_CI_all", scens, samplesizes, models)
saveRDS(cate_preds, paste0("live/results/continuous/all/cate_preds.RDS"))
print("cate done!")

BLP_folds <- extract_results("BLP_tests_all", scens, samplesizes, models)
saveRDS(BLP_folds, paste0("live/results/continuous/all/BLP_folds.RDS"))
print("BLP folds done!")

BLP_wholes <- extract_results("BLP_whole_all", scens, samplesizes, models)
saveRDS(BLP_wholes, paste0("live/results/continuous/all/BLP_wholes.RDS"))
print("BLP whole done!")

TEVIM_preds <- extract_results("te_vims_all", scens, samplesizes, models)
saveRDS(TEVIM_preds, paste0("live/results/continuous/all/te_vims.RDS"))
print("TEVIMS done!")

# draws takes a million years and makes the memory huge - leave for now
#draws_all <- extract_results("draws_all", scens, samplesizes, models[-1])