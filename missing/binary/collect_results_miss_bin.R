#####################
# title: collect up all the missing_binary results
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
scenarios <- c(1, 3, 8, 10)
sample_sizes <- c(250, 500, 1000)
model_names <- c("causal_forest", "dr_random_forest", "dr_superlearner")
result_types <- c("tau", "BLP_tests", "BLP_whole", "independence_tests", "independence_whole", "te_vims", "data", "truth")
miss_type <- c("prognostic", "predictive", "both")
miss_prop <- c(0.2, 0.4)
miss_method <- c("complete_cases", "mean_imputation", "multiple_imputation")

results_by_type <- setNames(vector("list", length(result_types)), result_types)

for (result_type in result_types) {
  results_by_type[[result_type]] <- list()
  for (scenario in scenarios) {
    scenario_name <- paste0("scenario_", scenario)
    results_by_type[[result_type]][[scenario_name]] <- list()
    
    for (n in sample_sizes) {
      sample_size <- paste0("size_", n)
      results_by_type[[result_type]][[scenario_name]][[sample_size]] <- list()
      
      for (type in miss_type) {
        results_by_type[[result_type]][[scenario_name]][[sample_size]][[type]] <- list()
        
        for (prop in miss_prop) {
          prop <- as.character(prop)
          results_by_type[[result_type]][[scenario_name]][[sample_size]][[type]][[prop]] <- list()
          
          for (method in miss_method) {
            results_by_type[[result_type]][[scenario_name]][[sample_size]][[type]][[prop]][[method]] <- list()
            
            res_dir <- paste0("live/results/missing/binary/", type, "/", prop, "/", method, "/scenario_", scenario, "/", n, "/all_methods")
            
            result_files <- list.files(res_dir, pattern = "res_sim_", full.names = TRUE)
            # To combine simulation runs for each model
            if (length(result_files) != 100) {
              print(paste0("sims missing for scenario ", scenario, " sample size ", n, " - skipping"))
              complete_sims <- list.files(res_dir, "res_sim")
              complete_nums <- gsub("res_sim_", "", complete_sims)
              complete_nums <- gsub(".RDS", "", complete_nums, ignore.case = T) %>% as.numeric()
              failed_sims <- setdiff(seq_len(100), complete_nums)
              print(paste0("number of missing sims: ", length(failed_sims)))
              print(failed_sims)
              next
            } 
            temp_model_res <- list()
            for (res_file in result_files) {
              sim_res <- readRDS(res_file)
              for (model in model_names) {
                value <- sim_res[[model]][[result_type]]
                if (!model %in% names(temp_model_res)) temp_model_res[[model]] <- list()
                temp_model_res[[model]] <- c(temp_model_res[[model]], list(value))
              }
            }
            # Optionally, simplify temp_model_res for each model (e.g., unlist or whatever fits your needs)
            results_by_type[[result_type]][[scenario_name]][[sample_size]][[type]][[prop]][[method]] <- temp_model_res
            
          }
        }
      }
    }
  }
}

# save the collected up results
saveRDS(results_by_type, paste0("live/results/new_format/competing_risk_all.RDS"))