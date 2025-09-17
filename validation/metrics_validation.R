#############################
# title: Compute metrics - CATE subgroup & variance validation
# updated: 21/08/2025
# author: Ellie Van Vogt
#############################

library(dplyr)
library(tidyr)

# Paths and parameters
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)
res_path <- file.path(path, "live/results")
results <- readRDS(file.path(res_path, "new_format/validation_all.RDS"))

scenarios    <- c(1, 2, 3)
sample_sizes <- c(250, 500, 1000)
interim_props <- c(0.25, 0.5, 0.75)
model_names  <- c("causal_forest", "dr_random_forest")

metrics <- list(
  subgroups = list(),
  variances = list()
)

#############
# SUBGROUPS METRICS
#############
for (scenario in scenarios) {
  scenario_name <- paste0("scenario_", scenario)
  metrics$subgroups[[scenario_name]] <- list()
  
  for (n in sample_sizes) {
    size_name <- paste0("size_", n)
    metrics$subgroups[[scenario_name]][[size_name]] <- list()
    
    for (iprop in interim_props) {
      iprop_name <- paste0("interim_", iprop)
      metrics$subgroups[[scenario_name]][[size_name]][[iprop_name]] <- list()
      
      # Get results object
      subgroup_list <- results$subgroups[[scenario_name]][[size_name]][[iprop_name]]
      if (is.null(subgroup_list)) next
      
      for (model in model_names) {
        model_list <- subgroup_list[[model]]
        if (is.null(model_list)) next
        model_metrics <- list()
        
        # Each element in the list is typically a named vector with top/bottom p-values
        for (sim_num in seq_along(model_list)) {
          vals <- model_list[[sim_num]]
          if (is.null(vals)) next
          sim_metrics <- data.frame(
            sim = sim_num,
            top_pval = as.numeric(vals[['top']]),
            bottom_pval = as.numeric(vals[['bottom']]),
            model = model,
            scenario = scenario,
            n = n,
            interim = iprop
          )
          model_metrics <- c(model_metrics, list(sim_metrics))
        }
        
        metrics$subgroups[[scenario_name]][[size_name]][[iprop_name]][[model]] <- bind_rows(model_metrics)
      }
    }
  }
}

#############
# VARIANCES METRICS
#############
for (scenario in scenarios) {
  scenario_name <- paste0("scenario_", scenario)
  metrics$variances[[scenario_name]] <- list()
  
  for (n in sample_sizes) {
    size_name <- paste0("size_", n)
    metrics$variances[[scenario_name]][[size_name]] <- list()
    
    for (iprop in interim_props) {
      iprop_name <- paste0("interim_", iprop)
      metrics$variances[[scenario_name]][[size_name]][[iprop_name]] <- list()
      
      variance_list <- results$variances[[scenario_name]][[size_name]][[iprop_name]]
      if (is.null(variance_list)) next
      
      for (model in model_names) {
        model_list <- variance_list[[model]]
        if (is.null(model_list)) next
        model_metrics <- list()
        
        # Each element in the list should be named vector c(vt1=..., vt2=...)
        for (sim_num in seq_along(model_list)) {
          vals <- model_list[[sim_num]]
          if (is.null(vals)) next
          sim_metrics <- data.frame(
            sim = sim_num,
            vt1 = as.numeric(vals['vt1']),
            vt2 = as.numeric(vals['vt2']),
            var_change = as.numeric(vals['vt2']) - as.numeric(vals['vt1']),
            model = model,
            scenario = scenario,
            n = n,
            interim = iprop
          )
          model_metrics <- c(model_metrics, list(sim_metrics))
        }
        
        metrics$variances[[scenario_name]][[size_name]][[iprop_name]][[model]] <- bind_rows(model_metrics)
      }
    }
  }
}

#############
# MERGE RESULTS TO DATAFRAMES
#############
# Subgroups
metrics_subgroups_df <- bind_rows(
  lapply(names(metrics$subgroups), function(scn) {
    lapply(names(metrics$subgroups[[scn]]), function(sz) {
      lapply(names(metrics$subgroups[[scn]][[sz]]), function(ip) {
        lapply(names(metrics$subgroups[[scn]][[sz]][[ip]]), function(mod) {
          metrics$subgroups[[scn]][[sz]][[ip]][[mod]]
        }) %>% bind_rows()
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) %>% bind_rows()
)

# Variances
metrics_variances_df <- bind_rows(
  lapply(names(metrics$variances), function(scn) {
    lapply(names(metrics$variances[[scn]]), function(sz) {
      lapply(names(metrics$variances[[scn]][[sz]]), function(ip) {
        lapply(names(metrics$variances[[scn]][[sz]][[ip]]), function(mod) {
          metrics$variances[[scn]][[sz]][[ip]][[mod]]
        }) %>% bind_rows()
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) %>% bind_rows()
)

#############
# SAVE OUTPUTS
#############
saveRDS(metrics,     file.path(res_path, "new_format/metrics_cate_subgroup_variance_nested.RDS"))
saveRDS(metrics_subgroups_df, file.path(res_path, "new_format/metrics_cate_subgroups_tidy.RDS"))
saveRDS(metrics_variances_df, file.path(res_path, "new_format/metrics_cate_variances_tidy.RDS"))
