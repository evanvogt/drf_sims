###############
# title: Compute performance metrics - competing risk methods
# updated: 21/08/2025
# author: Ellie Van Vogt
###############

library(dplyr)
library(tidyr)

# paths
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

res_path <- file.path(path, "live/results")

# arguments and parameters
scenarios     <- seq(1:7)
sample_sizes  <- c(250, 500, 1000)
model_names   <- c("causal_survival_forest_sub_dist", "causal_survival_forest_cause_spec")

# load collected results
results <- readRDS(file.path(res_path, "new_format/competing_risk_no_censoring.RDS"))

# Map each model to its corresponding truth column
model_truth_map <- list(
  causal_survival_forest_sub_dist = "tau_sdh",
  causal_survival_forest_cause_spec = "tau_cs1"    # use tau_cs2 if evaluating a different cause
)

metrics <- list()

for (scenario in scenarios) {
  scenario_name <- paste0("scenario_", scenario)
  metrics[[scenario_name]] <- list()
  
  for (n in sample_sizes) {
    size_name <- paste0("size_", n)
    metrics[[scenario_name]][[size_name]] <- list()
    
    # Get truth dataframe for this scenario & sample size
    truth_list <- results[["truth"]][[scenario_name]][[size_name]][["truth"]]
    
    if (is.null(truth_list)) {
      warning(paste("No truth_df found for scenario", scenario, "sample size", n, "- skipping"))
      next
    }
    
    for (model in model_names) {
      truth_col <- model_truth_map[[model]]
      
      temp_metrics <- list()
      for (sim_num in 1:1000) {
        # subset truth for this simulation
        truth_sim <- tryCatch(truth_list[[sim_num]], error = function(e) NULL)
        tau_est_obj <- tryCatch(results[["tau"]][[scenario_name]][[size_name]][[model]][[sim_num]], error = function(e) NULL)
        
        truth_val <- tryCatch(truth_sim[,truth_col], error = function(e) NULL)
        
        if (is.null(truth_sim) | is.null(tau_est_obj) | is.null(truth_val)) next
        
        tau_est <- tau_est_obj$tau
        
        bias <- mean(tau_est - truth_val, na.rm = T)
        corr <- cor(truth_val, tau_est, use = "pairwise.complete.obs")
        mse <- mean((tau_est - truth_val)^2)
        
        value <- data.frame(sim = sim_num, bias = bias, corr = corr, mse = mse,
                            model = model, scenario = scenario, n = n)
        temp_metrics <- c(temp_metrics, list(value))
      }
      
      # bind sims for model
      metrics[[scenario_name]][[size_name]][[model]] <- bind_rows(temp_metrics)
    }
  }
}

# Merge all results into one tidy dataframe
metrics_df <- bind_rows(
  lapply(names(metrics), function(scn) {
    lapply(names(metrics[[scn]]), function(sz) {
      lapply(names(metrics[[scn]][[sz]]), function(mod) {
        metrics[[scn]][[sz]][[mod]]
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) %>% bind_rows()
)

# Save outputs
saveRDS(metrics,     file.path(res_path, "new_format/metrics_competing_risks_no_cens_nested.RDS"))
saveRDS(metrics_df, file.path(res_path, "new_format/metrics_competing_risks_no_cens_tidy.RDS"))
