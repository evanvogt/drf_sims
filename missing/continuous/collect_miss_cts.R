#####################
# title: collect up all the missing_continuous results
####################
#libraries
library(dplyr)
library(tidyr)

# paths
base_path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"

# arguments and parameters
scenarios <- 3
sample_sizes <- 1000
model_names <- c("causal_forest", "dr_random_forest", "dr_superlearner")
result_types <- c("tau", "BLP_whole", "independence_whole", "data", "truth")

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
          prop_chr <- as.character(prop)
          results_by_type[[result_type]][[scenario_name]][[sample_size]][[type]][[prop_chr]] <- list()
          
          for (method in miss_method) {
            results_by_type[[result_type]][[scenario_name]][[sample_size]][[type]][[prop_chr]][[method]] <- list()
            
            res_dir <- file.path(base_path, "live/results/missing/continuous", 
                                 type, prop_chr, method, 
                                 paste0("scenario_", scenario), n, "all_methods")
            
            result_files <- list.files(res_dir, pattern = "^res_sim_\\d+\\.RDS$", full.names = TRUE)
            result_files <- result_files[order(as.numeric(gsub(".*res_sim_(\\d+)\\.RDS$", "\\1", 
                                                               result_files, ignore.case = TRUE)))]
            
            if (length(result_files) < 100) {
              warning(sprintf("Scenario %s, n=%s, type=%s, prop=%s, method=%s: only %d sims found (expected 100)", 
                              scenario, n, type, prop_chr, method, length(result_files)))
            }
            if (length(result_files) == 0) next
            
            temp_model_res <- list()
            for (res_file in result_files) {
              sim_res <- readRDS(res_file)
              
              if (result_type %in% c("data", "truth")) {
                sim_num <- gsub(".*res_sim_(\\d+)\\.RDS$", "\\1", res_file)
                value <- sim_res[[result_type]]
                if (!result_type %in% names(temp_model_res)) temp_model_res[[result_type]] <- list()
                temp_model_res[[result_type]] <- c(temp_model_res[[result_type]], list(value))
              } else {
                for (model in model_names) {
                  value <- sim_res[[model]][[result_type]]
                  if (!model %in% names(temp_model_res)) temp_model_res[[model]] <- list()
                  temp_model_res[[model]] <- c(temp_model_res[[model]], list(value))
                }
              }
            }
            results_by_type[[result_type]][[scenario_name]][[sample_size]][[type]][[prop_chr]][[method]] <- temp_model_res
          }
        }
      }
    }
  }
}

###############################
# Identify missing simulations (SCENARIO 3 only)
###############################
missing_array_ids <- c()
for (type in miss_type) {
  for (prop in miss_prop) {
    for (method in miss_method) {
      res_dir <- file.path(base_path, "live/results/missing/continuous", 
                           type, as.character(prop), method, 
                           "scenario_3", 1000, "all_methods")
      
      result_files <- list.files(res_dir, pattern = "^res_sim_\\d+\\.RDS$", full.names = TRUE)
      if (length(result_files) != 100) {
        complete_sims <- list.files(res_dir, "^res_sim_\\d+\\.RDS$")
        complete_nums <- gsub("res_sim_", "", complete_sims)
        complete_nums <- gsub("\\.RDS$", "", complete_nums, ignore.case = TRUE) %>% as.numeric()
        failed_sims <- setdiff(seq_len(100), complete_nums)
        
        # Index calculation for new array (matches new PBS script: 3*3*2*100 = 1800 jobs)
        # miss_method_idx: 0,1,2   miss_type_idx: 0,1,2   miss_prop_idx: 0,1 
        miss_type_idx <- match(type, miss_type) - 1
        miss_method_idx <- match(method, miss_method) - 1
        miss_prop_idx <- match(prop, miss_prop) - 1
        
        failed_array_indices <- (miss_method_idx * 600) + 
          (miss_type_idx * 200) +
          (miss_prop_idx * 100) +
          (failed_sims - 1) + 1
        missing_array_ids <- c(missing_array_ids, failed_array_indices)
      }
    }
  }
}

# save the list of missing array ids
failed_ids_path <- file.path(base_path, "live/scripts/drf_sims/missing/continuous/failed_ids.txt")
dir.create(dirname(failed_ids_path), recursive = TRUE, showWarnings = FALSE)
if (length(missing_array_ids) != 0) {
  cat(missing_array_ids, file = failed_ids_path, sep="\n")
}

# save the collected up results
results_outfile <- file.path(base_path, "live/results/new_format/missing_continuous_all.RDS")
dir.create(dirname(results_outfile), recursive = TRUE, showWarnings = FALSE)
saveRDS(results_by_type, results_outfile)

###############################
# Summary printout
###############################
cat("\n========== Summary of collected results ==========\n")

for (result_type in names(results_by_type)) {
  cat("\n--- Result type:", result_type, "---\n")
  for (scenario in names(results_by_type[[result_type]])) {
    for (sample_size in names(results_by_type[[result_type]][[scenario]])) {
      for (type in names(results_by_type[[result_type]][[scenario]][[sample_size]])) {
        for (prop in names(results_by_type[[result_type]][[scenario]][[sample_size]][[type]])) {
          for (method in names(results_by_type[[result_type]][[scenario]][[sample_size]][[type]][[prop]])) {
            
            temp <- results_by_type[[result_type]][[scenario]][[sample_size]][[type]][[prop]][[method]]
            
            if (length(temp) == 0) {
              cat(sprintf("Scenario %s, %s, %s, prop=%s, method=%s: 0 sims\n", 
                          scenario, sample_size, type, prop, method))
            } else {
              if (result_type %in% c("data", "truth")) {
                count <- length(temp[[result_type]])
                cat(sprintf("Scenario %s, %s, %s, prop=%s, method=%s: %d sims (data/truth)\n", 
                            scenario, sample_size, type, prop, method, count))
              } else {
                counts <- sapply(temp, length)
                counts_str <- paste(names(counts), counts, sep=": ", collapse=", ")
                cat(sprintf("Scenario %s, %s, %s, prop=%s, method=%s -> %s\n", 
                            scenario, sample_size, type, prop, method, counts_str))
              }
            }
          }
        }
      }
    }
  }
}
cat("\n=================================================\n")
