###############
# title: testing how much HPC resource is needed for one simulation run per model and scenario
# date started: 07/01/2025
# date finished:
# author: Ellie Van Vogt
###############
set.seed(1998)

# libraries ----
library(dplyr)
library(syrup)
library(hms)
library(doParallel)

# paths ----
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# parameters -----


samplesizes <- c(250, 500, 1000, 5000)
models <- c("CF", "DR_RF", "T_RF")# , "Logistic", "H_lasso") # not doing the last two because they are run by scenario anyways
scens <- paste0("scenario_", seq_along(1:10))
test_scenarios <- paste0("test_scenario_", seq_along(1:10))

# make some test datasets that are still lists to check the parameters for one simulation
for (scenario in scens) {
  for (n in samplesizes) {
    datasets <- readRDS(paste0(c("live/data/", scenario, "_", n, ".rds"), collapse = ""))
    
    test_dataset <- list(datasets[[1]])
    
    # now save as test_scenario_S_n
    saveRDS(test_dataset, paste0(c("live/data/test_", scenario, "_", n, ".rds"), collapse = ""))
  }
}

# also need to make some test directories that will deleted after the run
for (scenario in test_scenarios){
  dir.create(paste0(c("live", "results", scenario), collapse = "/"))
  
  for (s in samplesizes){
    dir.create(paste0(c("live", "results", scenario, s), collapse = "/"))

    for (m in models){
      dir.create(paste0(c("live", "results", scenario, s, m), collapse = "/"))
      
    }
  }
}
rm("datasets", "test_dataset", "m", "n", "s")
gc()


# make a function for profiling
profiling <- function(scenario) {
  samplesize_runs <- lapply(samplesizes, function(n) {
    args <- c(scenario, n, 1)
    assign("commandArgs", function(trailingOnly = TRUE) args, envir = globalenv())
    
    model_runs <- lapply(models, function(model) {
      model_scripts <- list.files(paste0("live/scripts/", model), ".R", full.names = T)
      model_scripts <- model_scripts[!grepl("test", model_scripts)]
      
      script_runs <- lapply(model_scripts, function(script) {
        usage <- syrup(
          source(script)
        )
        
        # extracting the values for execution of one dataset - this will be adapted in the script making stage for 1000 datasets
        time <- as_hms(max(usage$time, na.rm = T) - min(usage$time, na.rm = T))
        cpu <- max(usage$pct_cpu, na.rm = T)
        mem <- max(usage$rss, na.rm = T)
        
        script_name <- gsub(paste0("live/scripts/", model, "/"), "", script)
        script_name <- gsub(".R", "", script_name)
        scenario_name <- gsub("test_", "", scenario)

        return(list(script = script_name, scenario = scenario_name, n = n, time = time, cpu = cpu, mem = mem))
      }) %>% simplify2array()
      return(script_runs)
    }) %>% simplify2array()
    return(model_runs)
  }) %>% simplify2array()
  return(samplesize_runs)
}

#do the profiling in parallel

cl <- makeCluster(5, type = "PSOCK")
registerDoParallel(cl)

resource_data <- foreach(scenario = test_scenarios,
                         .packages = c("dplyr", "syrup", "hms")) %dopar% {
                           profiling(scenario)
                         }
stopCluster(cl)

#mclapply(test_scenarios, profiling, 1)#, mc.preschedule = F) %>% simplify2array()


# Restore the original commandArgs function from base R
assign("commandArgs", base::commandArgs, envir = globalenv())

#might be better to save stuff as an RDS to preserve all the data types?
saveRDS(resource_data, "live/data/resource_estimates.RDS")

# remove the test directories and test data files
for (scenario in test_scenarios) {
  unlink(paste0("live/results/", scenario), recursive = T)
}
test_dataset_list <- list.files("live/data", "test_scenario", full.names = T)
unlink(test_dataset_list)
