###############
# title: collect all the data
# date started: 20/02/25
# date finished:
# author: Ellie Van Vogt
###############

rm(list = ls(all = TRUE))

library(dplyr)
library(tidyr)

# paths
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# lists and categories
samplesizes <- c(100, 250, 500, 1000)
scens <- paste0("scenario_", seq_along(1:10))

data_list <- setNames(vector("list", length(scens)), scens)
for (scenario in scens) {
  data_list[[scenario]] <- setNames(vector("list", length(samplesizes)), as.character(samplesizes))
  
  for (n in samplesizes) {
    data <- readRDS(paste0("live/data/binary/", scenario, "_", n, ".RDS"))
    data_list[[scenario]][[as.character(n)]] <- data
  }
}

saveRDS(data_list, "live/data/binary/all_data.RDS")
