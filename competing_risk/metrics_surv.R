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
scenarios     <- c(1, 2, 3, 4)
sample_sizes  <- c(250, 500, 1000)
model_names   <- c("causal_survival_forest_sub_dist", "causal_survival_forest_cause_spec")

# load collected results
results <- readRDS(file.path(res_path, "new_format/competing_risk_all.RDS"))

# Map each model to its corresponding truth column
model_truth_map <- list(
  causal_survival_forest_sub_dist = "tau_subdist",
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
    truth_df <- results[["truth"]][[scenario_name]][[size_name]][["truth_df"]]
    
    if (is.null(truth_df)) {
      warning(paste("No truth_df found for scenario", scenario, "sample size", n, "- skipping"))
      next
    }
    
    for (model in model_names) {
      truth_col <- model_truth_map[[model]]
      
      if (!truth_col %in% colnames(truth_df)) {
        warning(paste("Truth column", truth_col, "not found for model", model,
                      "in scenario", scenario, "n", n))
        next
      }
      
      model_metrics <- list()
      
      for (sim_num in 1:100) {
        # subset truth for this simulation
        truth_sim <- truth_df %>% filter(sim == sim_num)
        
        if (nrow(truth_sim) == 0) {
          warning(paste("No truth rows for sim", sim_num,
                        "scenario", scenario, "n", n, "model", model))
          next
        }
        
        truth_vals <- truth_sim[[truth_col]]
        
        # retrieve tau estimates for this sim
        tau_obj <- results[["tau"]][[scenario_name]][[size_name]][[model]][[sim_num]]
        
        if (is.null(tau_obj)) {
          warning(paste("No estimates for sim", sim_num,
                        "scenario", scenario, "n", n, "model", model))
          next
        }
        
        tau_est <- tau_obj$tau
        
        # length check
        if (length(tau_est) != length(truth_vals)) {
          warning(paste("Length mismatch:",
                        "truth =", length(truth_vals),
                        "est =", length(tau_est),
                        "scenario", scenario, "n", n, "model", model, "sim", sim_num))
          next
        }
        
        # compute performance metrics
        bias <- mean(tau_est - truth_vals, na.rm = TRUE)
        
        # correlation only meaningful if tau varies (avoid polychoric NA)
        corr <- if (length(unique(truth_vals)) > 1) {
          suppressWarnings(cor(truth_vals, tau_est, use = "pairwise.complete.obs"))
        } else {
          NA_real_
        }
        
        mse  <- mean((tau_est - truth_vals)^2, na.rm = TRUE)
        
        sim_metrics <- data.frame(
          sim      = sim_num,
          bias     = bias,
          corr     = corr,
          mse      = mse,
          estimand = truth_col,
          model    = model,
          scenario = scenario,
          n        = n
        )
        
        model_metrics <- c(model_metrics, list(sim_metrics))
      }
      
      # bind sims for model
      metrics[[scenario_name]][[size_name]][[model]] <- bind_rows(model_metrics)
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
saveRDS(metrics,     file.path(res_path, "new_format/metrics_competing_risks_nested.RDS"))
saveRDS(metrics_df, file.path(res_path, "new_format/metrics_competing_risks_tidy.RDS"))
