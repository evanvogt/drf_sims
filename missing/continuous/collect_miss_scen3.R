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
