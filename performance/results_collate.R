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
source("live/scripts/functions/extract_results.R")

# headings
# i thjikn results up to 1000 will be available for all the models fingers crossed
samplesizes <- c(250, 500, 1000, 5000)
models <- c("CF", "DR_RF", "DR_oracle", "T_RF")#, "Logistic", "H_lasso")
scens <- paste0("scenario_", seq_along(1:10))

# combine results by type I guess so that it's easier to do the performance metrics

res_types <- c("tau", "BLP_tests", "BLP_whole", "te_vims", )

cate_preds <- extract_results("tau", scens, samplesizes, models)
BLP_folds <- extract_results("BLP_tests", scens, samplesizes, models)
BLP_wholes <- extract_results("BLP_whole", scens, samplesizes, models)
TEVIM_preds <- extract_results("te_vims", scens, samplesizes, models)
draws_all <- extract_results("draws", scens, samplesizes, models[-1])


# save this mega list of results
res_files <- c("cate_preds", "BLP_folds", "BLP_wholes", "TEVIM_preds", "draws_all")

for (res in res_files) {
  saveRDS(as.name(res), paste0("live/results/all/", res, ".RDS"))
}
