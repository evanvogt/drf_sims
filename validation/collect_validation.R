#####################
# title: collect up all validation results (simplified tidy)
# date started: 21/08/2025
# author: Ellie Van Vogt
#####################

library(dplyr)
library(tidyr)

# paths
base_path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"

# arguments and parameters
scenarios <- 3        
sample_sizes <- 1000
interim_props <- c(0.25, 0.5, 0.75) 
model_names <- c("causal_forest", "dr_random_forest")
validation_types <- c("subgroups", "variances", "var_imps")

# Initialize list to hold tidy results
tidy_results <- setNames(vector("list", length(validation_types)), validation_types)

#####################
# Main loop: Collect into tidy dataframes
#####################
for (val_type in validation_types) {
  res_rows <- list()
  
  for (scenario in scenarios) {
    scenario_name <- paste0("scenario_", scenario)
    
    for (n in sample_sizes) {
      size_name <- paste0("size_", n)
      
      for (iprop in interim_props) {
        iprop_name <- paste0("interim_", iprop)
        
        res_dir <- file.path(
          base_path,
          "live/results/validation",
          paste0("scenario_", scenario),
          as.character(n),
          as.character(iprop)
        )
        
        result_files <- list.files(res_dir, pattern = "res_sim_", full.names = TRUE)
        
        if (length(result_files) == 0) next
        
        for (sim_num in seq_along(result_files)) {
          res_file <- result_files[sim_num]
          sim_res <- readRDS(res_file)
          value <- sim_res$validations[[val_type]]
          
          if (val_type %in% c("subgroups", "variances", "var_imps")) {
            for (model in model_names) {
              model_val <- value[[model]]
              # Subgroups: expect named vector with 'top' and 'bottom'
              if (val_type == "subgroups" && !is.null(model_val)) {
                res_rows[[length(res_rows)+1]] <- data.frame(
                  scenario = scenario,
                  sample_size = n,
                  interim_prop = iprop,
                  sim = sim_num,
                  model = model,
                  top_pval = as.numeric(model_val["top"]),
                  bottom_pval = as.numeric(model_val["bottom"])
                )
              }
              # Variances: expect named vector with vt1 and vt2
              if (val_type == "variances" && !is.null(model_val)) {
                res_rows[[length(res_rows)+1]] <- data.frame(
                  scenario = scenario,
                  sample_size = n,
                  interim_prop = iprop,
                  sim = sim_num,
                  model = model,
                  vt1 = as.numeric(model_val["vt1"]),
                  vt2 = as.numeric(model_val["vt2"]),
                  var_change = as.numeric(model_val["vt2"]) - as.numeric(model_val["vt1"])
                )
              }
              # Var imps: repeat every 4 (names, vi1, vi2, diff) – usually only 1 per sim
              if (val_type == "var_imps" && !is.null(model_val)) {
                n_entries <- length(model_val) / 4
                for (k in seq_len(n_entries)) {
                  idx <- (k-1)*4
                  varnames <- model_val[[idx+1]]
                  vi1 <- model_val[[idx+2]]
                  vi2 <- model_val[[idx+3]]
                  diff <- model_val[[idx+4]]
                  df <- data.frame(
                    scenario = scenario,
                    sample_size = n,
                    interim_prop = iprop,
                    sim = sim_num,
                    model = model,
                    variable = varnames,
                    vi1 = vi1,
                    vi2 = vi2,
                    diff = diff
                  )
                  res_rows[[length(res_rows)+1]] <- df
                }
              }
            } # end model loop
          }
        } # end sim loop
      }
    }
  }
  # Combine all rows to a single tidy dataframe per type
  tidy_results[[val_type]] <- bind_rows(res_rows)
}

# Save results
out_dir <- file.path(base_path, "live/results/new_format")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(tidy_results, file.path(out_dir, "validation_all_tidy.RDS"))

#####################
# Summary printout
#####################
cat("\n========== Summary of collected validation results ==========\n")
for (val_type in names(tidy_results)) {
  nrows <- nrow(tidy_results[[val_type]])
  cat(sprintf("Validation type %s: %d rows collected\n", val_type, nrows))
}
cat("\n=================================================\n")
