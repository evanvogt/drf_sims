#####################
# title: collect up all the binary results
# date started: 18/03/2025
# date finished:
# author: Ellie Van Vogt
####################
#libraries
library(dplyr)
library(tidyr)

# paths
base_path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"

# arguments and parameters
scenarios <- c(1:10)
sample_sizes <- c(100, 250, 500, 1000)
model_names <- c("causal_forest", "dr_random_forest", "dr_superlearner", "dr_oracle", "dr_semi_oracle")
result_types <- c("tau", "BLP_whole", "independence_whole", "data", "truth")

results_by_type <- setNames(vector("list", length(result_types)), result_types)

for (result_type in result_types) {
  results_by_type[[result_type]] <- list()
  
  for (scenario in scenarios) {
    scenario_name <- paste0("scenario_", scenario)
    results_by_type[[result_type]][[scenario_name]] <- list()
    
    for (n in sample_sizes) {
      sample_size <- paste0("size_", n)
      results_by_type[[result_type]][[scenario_name]][[sample_size]] <- list()
      
      res_dir <- file.path(base_path, "live/results/binary", 
                           paste0("scenario_", scenario), as.character(n), "all_methods")
      
      result_files <- list.files(res_dir, pattern = "res_sim_", full.names = TRUE)
      
      if (length(result_files) < 1000) {
        warning(sprintf("Scenario %s, size %s: only %d sims found (expected 1000)", 
                        scenario, n, length(result_files)))
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
      results_by_type[[result_type]][[scenario_name]][[sample_size]] <- temp_model_res
    }
  }
}

# Get the IDs for the failed simulations, organized by sample size:
missing_ids_by_sample_size <- list()

for (n in sample_sizes) {
  missing_ids_by_sample_size[[as.character(n)]] <- c()
}

for (scenario in scenarios) {
  for (n in sample_sizes) {
    res_dir <- file.path(base_path, "live/results/binary", 
                         paste0("scenario_", scenario), as.character(n), "all_methods")
    
    result_files <- list.files(res_dir, pattern = "res_sim_", full.names = TRUE)
    if (length(result_files) != 1000) {
      complete_sims <- list.files(res_dir, "res_sim")
      complete_nums <- gsub("res_sim_", "", complete_sims)
      complete_nums <- gsub(".RDS", "", complete_nums, ignore.case = TRUE) %>% as.numeric()
      failed_sims <- setdiff(seq_len(100), complete_nums)
      sample_idx <- match(n, sample_sizes) - 1
      scenario_idx <- match(scenario, scenarios) - 1
      failed_array_indices <- (sample_idx * 400) + (scenario_idx * 100) + (failed_sims - 1) + 1
      missing_ids_by_sample_size[[as.character(n)]] <- c(missing_ids_by_sample_size[[as.character(n)]], failed_array_indices)
    }
  }
}

# Sort each list
for (n in sample_sizes) {
  missing_ids_by_sample_size[[as.character(n)]] <- sort(missing_ids_by_sample_size[[as.character(n)]])
}
# Save the lists of missing array ids separately by sample size
failed_ids_dir <- file.path(base_path, "live/scripts/drf_sims/binary/all_methods/jobscripts")
dir.create(failed_ids_dir, recursive = TRUE, showWarnings = FALSE)

for (n in sample_sizes) {
  missing_ids <- missing_ids_by_sample_size[[as.character(n)]]
  if (length(missing_ids) != 0) {
    failed_ids_path <- file.path(failed_ids_dir, paste0("failed_ids_", n, ".txt"))
    cat(missing_ids, file = failed_ids_path, sep = "\n")
  }
}

# save the collected up results
results_outfile <- file.path(base_path, "live/results/new_format/binary_all.RDS")
dir.create(dirname(results_outfile), recursive = TRUE, showWarnings = FALSE)
saveRDS(results_by_type, results_outfile)

#####################
# Summary printout
#####################
cat("\n========== Summary of collected results ==========\n")

for (result_type in names(results_by_type)) {
  cat("\n--- Result type:", result_type, "---\n")
  for (scenario in names(results_by_type[[result_type]])) {
    for (sample_size in names(results_by_type[[result_type]][[scenario]])) {
      temp <- results_by_type[[result_type]][[scenario]][[sample_size]]
      
      if (length(temp) == 0) {
        cat(sprintf("Scenario %s, %s: 0 simulations collected\n", scenario, sample_size))
      } else {
        if (result_type %in% c("data", "truth")) {
          count <- length(temp[[result_type]])
          cat(sprintf("Scenario %s, %s: %d simulations (data/truth)\n", 
                      scenario, sample_size, count))
        } else {
          counts <- sapply(temp, length)
          counts_str <- paste(names(counts), counts, sep=": ", collapse=", ")
          cat(sprintf("Scenario %s, %s -> %s\n", scenario, sample_size, counts_str))
        }
      }
    }
  }
}

cat("\n=================================================\n")
