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

# headings
samplesizes <- c(250, 500, 1000, 5000)
models <- c("CF", "DR_RF", "T_RF", "Logistic", "H_lasso")
scens <- paste0("scenario_", seq_along(1:10))


# combine all the results across all of the scenarios

# structure: scenario -> sample size -> model -> results (named by filename)
sim_results <- setNames(vector("list", length(scens)), scens)

for (scenario in scens) {
  sim_results[[scenario]] <- setNames(vector("list", length(samplesizes)), as.character(samplesizes))
  
  for (n in samplesizes) {
    sim_results[[scenario]][[as.character(n)]] <- setNames(vector("list", length(models)), models)
    
    for (model in models) {
      dir_path <- paste0("live/results/", scenario, "/", n, "/", model)
      
      if (dir.exists(dir_path)) {
        files <- list.files(dir_path, full.names = TRUE)
        
        if (length(files) > 0) {  # check there are results to put in
          filenames <- gsub(".RDS", "", files)
          res_list <- setNames(lapply(files, readRDS), basename(filenames))  # Named list: filename -> result
          
          sim_results[[scenario]][[as.character(n)]][[model]] <- res_list
        }
      }
    }
  }
}

# save this mega list of results
saveRDS(sim_results, "live/results/all/sim_results_all.RDS")

# maybe it's a good idea to collect all the datasets into one list as well?


data_path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/data"
data_all <- setNames(vector("list", length(scens)), scens)
data_list <- list.files(data_path)

for (scenario in scens) {
  data_all[[scenario]] <- setNames(vector("list", length(samplesizes)), as.character(samplesizes))
  
  for (n in samplesizes) {
    # Construct expected filename: scenario_X_Y (e.g., scenario_1_250)
    file_name <- data_list[grepl(paste0(scenario, "_", n, ".rds"), data_list, ignore.case = T)]
    file_path <- file.path(data_path, file_name)
    
    if (file.exists(file_path)) {
      data_all[[scenario]][[as.character(n)]] <- readRDS(file_path)
    }
  }
}

saveRDS(data_all, paste0(data_path, "/all_data.RDS"))
