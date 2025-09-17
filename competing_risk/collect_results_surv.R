#####################
# title: Collect competing_risk results (fixed, with proper truth handling + summary)
# updated: 21/08/2025
# author: Ellie Van Vogt
####################

library(dplyr)
library(tidyr)

# Set working directory
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# Parameters
scenarios     <- c(1, 2, 3, 4)
sample_sizes  <- c(250, 500, 1000)
model_names   <- c("causal_survival_forest_sub_dist", "causal_survival_forest_cause_spec")
result_types  <- c("tau", "data", "truth")

results_by_type <- setNames(vector("list", length(result_types)), result_types)

for (result_type in result_types) {
  results_by_type[[result_type]] <- list()
  
  for (scenario in scenarios) {
    scenario_name <- paste0("scenario_", scenario)
    results_by_type[[result_type]][[scenario_name]] <- list()
    
    for (n in sample_sizes) {
      size_name <- paste0("size_", n)
      results_by_type[[result_type]][[scenario_name]][[size_name]] <- list()
      
      res_dir <- file.path("live/results/competing_risk", scenario_name, n, "all_methods")
      result_files <- list.files(res_dir, pattern = "res_sim_", full.names = TRUE)
      
      # Only proceed if all simulation files are present
      if (length(result_files) != 100) {
        message("sims missing for ", scenario_name, " sample size ", n, " - skipping")
        next
      }
      
      # Order files by simulation number for consistent binding
      sim_nums <- as.integer(gsub(".*res_sim_(\\d+)\\.RDS$", "\\1", basename(result_files)))
      file_order <- order(sim_nums)
      result_files <- result_files[file_order]
      sim_nums <- sim_nums[file_order]
      
      temp_res <- list()
      
      for (i in seq_along(result_files)) {
        res_file <- result_files[i]
        sim_num  <- sim_nums[i]
        sim_res  <- readRDS(res_file)
        
        if (result_type == "truth") {
          # Collect the full truth data frame for each simulation
          value <- sim_res[["truth"]]
          if (is.null(value)) next
          value$sim <- sim_num
          
          if (is.null(temp_res[["truth_df"]])) {
            temp_res[["truth_df"]] <- value
          } else {
            temp_res[["truth_df"]] <- bind_rows(temp_res[["truth_df"]], value)
          }
          
        } else if (result_type == "data") {
          # Collect data as a list (by simulation number)
          value <- sim_res[["data"]]
          if (is.null(temp_res[["data_list"]])) temp_res[["data_list"]] <- vector("list", 100)
          temp_res[["data_list"]][[sim_num]] <- value
          
        } else if (result_type == "tau") {
          # Collect tau for each model and simulation
          for (model in model_names) {
            if (is.null(temp_res[[model]])) temp_res[[model]] <- vector("list", 100)
            temp_res[[model]][[sim_num]] <- sim_res[[model]][["tau"]]
          }
        }
      }
      
      results_by_type[[result_type]][[scenario_name]][[size_name]] <- temp_res
    }
  }
}

## Identify missing simulations (for array job tracking)
missing_array_ids <- c()
for (scenario in scenarios) {
  for (n in sample_sizes) {
    res_dir <- file.path("live/results/competing_risk", paste0("scenario_", scenario), n, "all_methods")
    result_files <- list.files(res_dir, pattern = "res_sim_", full.names = TRUE)
    
    if (length(result_files) != 100) {
      complete_sims <- gsub("res_sim_(\\d+)\\.RDS", "\\1", list.files(res_dir, pattern = "res_sim_"))
      complete_nums <- as.integer(complete_sims)
      failed_sims <- setdiff(seq_len(100), complete_nums)
      
      sample_idx   <- match(n, sample_sizes) - 1
      scenario_idx <- match(scenario, scenarios) - 1
      failed_array_indices <- (sample_idx * 400) + (scenario_idx * 100) + (failed_sims - 1) + 1
      
      missing_array_ids <- c(missing_array_ids, failed_array_indices)
    }
  }
}

if (length(missing_array_ids) > 0) {
  cat(missing_array_ids,
      file = "live/scripts/drf_sims/competing_risk/failed_ids.txt",
      sep = "\n")
}

# Save the collected results for analysis
saveRDS(results_by_type, "live/results/new_format/competing_risk_all.RDS")

#####################
# Summary printout
#####################
cat("\n========== Summary of collected results ==========\n")

for (result_type in names(results_by_type)) {
  cat("\n--- Result type:", result_type, "---\n")
  for (scenario in names(results_by_type[[result_type]])) {
    for (size_name in names(results_by_type[[result_type]][[scenario]])) {
      temp <- results_by_type[[result_type]][[scenario]][[size_name]]
      
      if (length(temp) == 0) {
        cat(sprintf("Scenario %s, %s: 0 simulations collected\n", scenario, size_name))
      } else {
        if (result_type == "truth") {
          # truth_df may be missing if not collected
          if (!is.null(temp[["truth_df"]])) {
            count <- length(unique(temp[["truth_df"]]$sim))
            cat(sprintf("Scenario %s, %s: %d simulations (truth_df)\n",
                        scenario, size_name, count))
          } else {
            cat(sprintf("Scenario %s, %s: 0 simulations (truth_df missing)\n",
                        scenario, size_name))
          }
        } else if (result_type == "data") {
          # data_list usually has 100 entries, some may be NULL
          if (!is.null(temp[["data_list"]])) {
            count <- sum(sapply(temp[["data_list"]], Negate(is.null)))
            cat(sprintf("Scenario %s, %s: %d simulations (data_list)\n",
                        scenario, size_name, count))
          } else {
            cat(sprintf("Scenario %s, %s: 0 simulations (data_list missing)\n",
                        scenario, size_name))
          }
        } else if (result_type == "tau") {
          # For tau, print out counts per model
          model_counts <- sapply(model_names, function(mod) {
            if (!is.null(temp[[mod]])) sum(sapply(temp[[mod]], Negate(is.null))) else 0
          })
          counts_str <- paste(model_names, model_counts, sep=": ", collapse=", ")
          cat(sprintf("Scenario %s, %s tau counts -> %s\n", scenario, size_name, counts_str))
        }
      }
    }
  }
}

if (length(missing_array_ids) > 0) {
  cat(sprintf("\n\nMissing simulation array IDs (%d total):\n", length(missing_array_ids)))
  cat(paste(missing_array_ids, collapse=", "), "\n")
}

cat("\n=================================================\n")
