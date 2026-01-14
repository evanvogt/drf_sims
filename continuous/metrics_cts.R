###############
# title: Compute performance metrics - continuous outcome
# updated 20/08/2025
###############

library(dplyr)
library(tidyr)

# paths
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

res_path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/results/"

# Ensure output directory exists
dir.create(file.path(res_path, "new_format"), showWarnings = FALSE, recursive = TRUE)

# arguments and parameters
scenarios <- c(1:10)
sample_sizes <- c(100, 250, 500, 1000)
model_names <- c("causal_forest", "dr_random_forest", "dr_superlearner", "dr_oracle", "dr_semi_oracle")

results <- readRDS(paste0(res_path, "new_format/continuous_all.RDS"))

threshold <- 0.1

metrics <- list()
for (scenario in scenarios) {
  scenario_name <- paste0("scenario_", scenario)
  metrics[[scenario_name]] <- list()
  
  for (n in sample_sizes) {
    sample_size <- paste0("size_", n)
    metrics[[scenario_name]][[sample_size]] <- list()
    
    # get truth data frame for this scenario/size
    truth <- results[["truth"]][[scenario_name]][[sample_size]][["truth"]]
    
    for (model in model_names) {
      metrics[[scenario_name]][[sample_size]][[model]] <- list()
      
      temp_metrics <- list()
      for (sim_num in 1:1000) {
        # Check if all necessary data exist for this simulation
        # If missing, skip this simulation
        
        truth_sim <- tryCatch(truth[[sim_num]], error = function(e) NULL)
        tau_est_obj <- tryCatch(results[["tau"]][[scenario_name]][[sample_size]][[model]][[sim_num]], error = function(e) NULL)
        BLP_obj <- tryCatch(results[["BLP_whole"]][[scenario_name]][[sample_size]][[model]][[sim_num]], error = function(e) NULL)
        indep_obj <- tryCatch(results[["independence_whole"]][[scenario_name]][[sample_size]][[model]][[sim_num]], error = function(e) NULL)
        
        if (is.null(truth_sim) | is.null(tau_est_obj) | is.null(BLP_obj) | is.null(indep_obj)) next
        
        truth_val <- truth_sim$tau
        tau_est <- tau_est_obj  # estimated tau vector
        BLP <- BLP_obj
        indep <- indep_obj
        
        BLP_p <- tryCatch({
          BLP[4,2]
        }, error = function(e) { 1 })
        indep_p <- indep$p_value
        
        indep_failed <- as.numeric(indep$method == "independence_test_failed")
        
        bias <- mean(tau_est - truth_val, na.rm = TRUE)
        corr <- ifelse(scenario != 1, cor(truth_val, tau_est, use = "pairwise.complete.obs"), 0)
        mse  <- mean((tau_est - truth_val)^2, na.rm = TRUE)
        
        if (scenario == 1) {
          BLP_correct <- ifelse(BLP_p > threshold, 1, 0)
          indep_correct <- ifelse(indep_p > threshold, 1, 0)
        } else {
          BLP_correct <- ifelse(BLP_p < threshold, 1, 0)
          indep_correct <- ifelse(indep_p < threshold, 1, 0)
        }
        
        value <- data.frame(sim = sim_num, bias = bias, corr = corr, mse = mse,
                            BLP_p = as.numeric(BLP_p), BLP_correct = BLP_correct,
                            indep_p = as.numeric(indep_p), indep_correct = indep_correct,
                            indep_failed = indep_failed, model = model, scenario = scenario, n = n)
        temp_metrics <- c(temp_metrics, list(value))
      }
      metrics[[scenario_name]][[sample_size]][[model]] <- bind_rows(temp_metrics)
    }
  }
}

# Combine into one tidy data frame for easier plotting/analysis
metrics_df <- bind_rows(
  lapply(names(metrics), function(scn) {
    lapply(names(metrics[[scn]]), function(sz) {
      lapply(names(metrics[[scn]][[sz]]), function(mod) {
        metrics[[scn]][[sz]][[mod]]
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) %>% bind_rows()
)

saveRDS(metrics, paste0(res_path, "new_format/metrics_continuous_nested.RDS"))
saveRDS(metrics_df, paste0(res_path, "new_format/metrics_continuous_tidy.RDS"))
