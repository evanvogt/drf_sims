##########
# Title: collect up the results - missing cts outcome
##########

# libraries
library(here)
library(dplyr)

# paths
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live"

# parameters
scenarios <- paste0("scenario_", seq(1:5))
sample_sizes <- as.character(500)
types <- c("prognostic", "predictive", "both")
props <- as.character(0.3)
mechanisms <- c("MAR", "AUX")
methods <- c("complete_cases", "mean_imputation", "missforest", "regression", "missing_indicator", "IPW", "none")


results <- list()
failed <- data.frame(matrix(nrow = 0, ncol = 7))
colnames(failed) <- c("scenario", "n", "prop", "type", "mechanism", "method", "run")

for (scenario in scenarios) {
  results[[scenario]] <- list()
  
  for (n in sample_sizes) {
    results[[scenario]][[as.character(n)]] <- list()
    
    for (type in types) {
      if (scenario == "scenario_1" & type != "prognostic") next
      
      results[[scenario]][[n]][[type]] <- list()
      
      for (prop in props) {
        results[[scenario]][[n]][[type]][[prop]] <- list()
        
        for (mechanism in mechanisms) {
          results[[scenario]][[n]][[type]][[prop]][[mechanism]] <- list()
          
          for (method in methods) {
            results[[scenario]][[n]][[type]][[prop]][[mechanism]][[method]] <- list()
            
            folder <- file.path(path, "results/missing/continuous", scenario, n, type, prop, mechanism, method)
            
            result_files <- list.files(folder, pattern = "^res_sim_\\d+\\.RDS$", full.names = TRUE)

            
            if (length(result_files) == 0) next
            
            temp <- list()
            for (res_file in result_files) {
              sim_res <- readRDS(res_file)
              sim_num <- gsub(".*res_sim_(\\d+)\\.RDS$", "\\1", res_file)
              temp[[sim_num]] <- sim_res
            }
            results[[scenario]][[n]][[type]][[prop]][[mechanism]][[method]] <- temp
            
            # are there sims missing?
            if (length(result_files) < 100) {
              complete_runs <- gsub(".*res_sim_(\\d+)\\.RDS$", "\\1", result_files) %>% as.numeric()
              failed_runs <- setdiff(seq_len(100), complete_runs)
              failed_params <- data.frame(scenario = gsub("^scenario_", "", scenario),
                                          n = as.numeric(n),
                                          prop = as.numeric(prop),
                                          type = type,
                                          mechanism = mechanism,
                                          method = method,
                                          run = failed_runs)
              failed <- rbind(failed, failed_params)
            }
          }
        }
      }
    }
  }
}

# save the big results file:
res_output <- file.path(path, "results/new_format/missing_continuous_all.RDS")
dir.create(dirname(res_output), recursive = T, showWarnings = F)
saveRDS(results, res_output)
print("Collection complete!")

# find array numbers for the failed simulations
if (nrow(failed) > 0) {
  params <- expand.grid(scenario = c(1:5),
                        n = c(500),
                        prop = c(0.3),
                        type = c("prognostic", "predictive", "both"),
                        mechanism = c("MAR", "AUX"),
                        method = c("complete_cases", "mean_imputation", "missforest", "regression", "missing_indicator", "IPW", "none"),
                        run = c(1:100))
  params <- params %>%
    filter(!(scenario == 1 & type != "prognostic"))
  
  failed_idx <- which(interaction(params) %in% interaction(failed))
  # save failed ids for rerunning
  failed_file <- file.path(path, "scripts/drf_sims/missing/continuous/jobscripts/failed_ids.txt")
  cat(failed_idx, file = failed_file, sep = "\n")
  print(paste0("failed runs found (", nrow(failed), ") saved to jobscripts folder"))
}