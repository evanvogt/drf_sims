#####################
# title: collecting results in an array job to make things quicker
# date started: 18/03/2025
# date finished:
# author: Ellie Van Vogt
####################

#libraries
library(dplyr)
library(tidyr)

# paths
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# arguments and parameters
args <- commandArgs(trailingOnly = T)
scenario <- as.character(args[1])
n <- as.numeric(args[2])
model <- as.character(args[3])

# path to results
res_dir <- paste0("live/results/", scenario, "/", n, "/", model)

# all the simulation runs
all_files <- list.files(res_dir, "res_sim", full.names = T)

if (length(all_files) == 1000) {
  res_sims_list <- lapply(all_files, readRDS)
  
  # all the taus and CIs
  tau_list <- lapply(res_sims_list, `[[`, "tau")
  saveRDS(tau_list, paste0(res_dir, "/tau_CI_all.RDS"))
  
  # BLP tests on each fold
  BLP_folds <- lapply(res_sims_list, `[[`, "BLP_tests")
  saveRDS(BLP_folds, paste0(res_dir, "/BLP_tests_all.RDS"))
  
  # BLP on the whole dataset
  BLP_wholes <- lapply(res_sims_list, `[[`, "BLP_whole")
  saveRDS(BLP_wholes, paste0(res_dir, "/BLP_whole_all.RDS"))
  
  # TEVIMS 
  te_vims_list <- lapply(res_sims_list, `[[`, "te_vims") 
  saveRDS(te_vims_list, paste0(res_dir, "/te_vims_all.RDS"))
  
  if (model != "CF") {
    draws_list <- lapply(res_sims_list, `[[`, "draws")
    saveRDS(draws_list, paste0(res_dir, "/draws_all.RDS"))
  }
  
  print(paste0(scenario, "_", n, " ", model, " has all the sims and results have been collated"))
  
} else {
  complete_sims <- list.files(res_dir, "res_sim")
  complete_nums <- gsub("res_sim_", "", complete_sims)
  complete_nums <- gsub(".RDS", "", complete_nums, ignore.case = T) %>% as.numeric()
  failed_sims <- setdiff(seq_len(1000), complete_nums)
  
  scen_num <- substr(scenario, nchar(scenario), nchar(scenario)) %>% as.numeric()
  array_nums <- 1000*(scen_num - 1) + failed_sims %>% sort()
  
  writeLines(as.character(array_nums), paste0("live/scripts/", model, "/jobscripts/failed_", n, "_", scenario, ".txt"))
  
  print(paste0(scenario, "_", n, " ", model, " has ", length(failed_sims), " missing sims, missing sims list saved to jobscripts folder"))
}

