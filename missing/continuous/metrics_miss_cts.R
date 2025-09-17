###############
# title: Compute performance metrics - continuous outcome with missingness
# updated 21/08/2025
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
scenarios <- 1:3 # corrected: should be a vector if you want >1 scenario
sample_sizes <- c(1000) # corrected: vector for sample sizes (add more if needed)
model_names <- c("causal_forest", "dr_random_forest", "dr_superlearner")
miss_type <- c("prognostic", "predictive", "both")
miss_prop <- c(0.2, 0.4)
miss_method <- c("complete_cases", "mean_imputation", "multiple_imputation")

results <- readRDS(paste0(res_path, "new_format/missing_continuous_all.RDS"))

threshold <- 0.1

metrics <- list()
for (scenario in scenarios) {
  scenario_name <- paste0("scenario_", scenario)
  metrics[[scenario_name]] <- list()
  
  for (n in sample_sizes) {
    sample_size <- paste0("size_", n)
    metrics[[scenario_name]][[sample_size]] <- list()
    
    for (type in miss_type) {
      metrics[[scenario_name]][[sample_size]][[type]] <- list()
      
      for (prop in miss_prop) {
        prop_name <- as.character(prop)
        metrics[[scenario_name]][[sample_size]][[type]][[prop_name]] <- list()
        
        for (method in miss_method) {
          metrics[[scenario_name]][[sample_size]][[type]][[prop_name]][[method]] <- list()
          
          # Checks in case the structure is missing
          if (!is.null(results[["truth"]][[scenario_name]][[sample_size]][[type]][[prop_name]][[method]][["truth"]])) {
            truth <- results[["truth"]][[scenario_name]][[sample_size]][[type]][[prop_name]][[method]][["truth"]]
          } else {
            next # skip this combination if missing
          }
          
          for (model in model_names) {
            temp_metrics <- list()
            for (sim_num in 1:100) {
              # Defensive: check if expected elements are present
              truth_sim <- tryCatch(truth[[sim_num]], error = function(e) NULL)
              if (is.null(truth_sim)) next
              
              truth_val <- truth_sim$tau
              # model estimates
              tau_est <- tryCatch(
                results[["tau"]][[scenario_name]][[sample_size]][[type]][[prop_name]][[method]][[model]][[sim_num]],
                error = function(e) NA
              )
              # Defensive: tau_est could be NULL or missing
              if (is.null(tau_est) || all(is.na(tau_est))) next
              # tau_est <- tau_est$tau # uncomment if needed
              
              # HTE test results
              BLP <- tryCatch(
                results[["BLP_whole"]][[scenario_name]][[sample_size]][[type]][[prop_name]][[method]][[model]][[sim_num]],
                error = function(e) NA
              )
              indep <- tryCatch(
                results[["independence_whole"]][[scenario_name]][[sample_size]][[type]][[prop_name]][[method]][[model]][[sim_num]],
                error = function(e) NA
              )
              
              BLP_p <- tryCatch({
                BLP[4,2]
              }, error = function(e) {
                1
              }) 
              indep_p <- tryCatch({
                indep$p_value
              }, error = function(e) {
                1
              })
              
              # metrics
              bias <- mean(tau_est - truth_val, na.rm = TRUE)
              corr <- ifelse(scenario != 1, cor(truth_val, tau_est, use = "pairwise.complete.obs"), 0)
              mse  <- mean((tau_est - truth_val)^2, na.rm = TRUE)
              
              if (scenario == 1) {
                BLP_correct <- ifelse(BLP_p > threshold, 1, 0)
                indep_correct <- ifelse(indep_p > threshold, 1, 0)
              } else {
                BLP_correct <- ifelse(BLP_p > threshold, 0, 1)
                indep_correct <- ifelse(indep_p > threshold, 0, 1)
              }
              
              value <- data.frame(
                sim = sim_num,
                bias = bias,
                corr = corr,
                mse = mse,
                BLP_correct = BLP_correct,
                indep_correct = indep_correct,
                model = model,
                scenario = scenario,
                n = n,
                miss_type = type,
                miss_prop = prop,
                miss_method = method,
                stringsAsFactors = FALSE
              )
              temp_metrics[[length(temp_metrics) + 1]] <- value
            }
            # Only bind if we actually got results
            if (length(temp_metrics) > 0) {
              metrics[[scenario_name]][[sample_size]][[type]][[prop_name]][[method]][[model]] <- bind_rows(temp_metrics)
            } else {
              metrics[[scenario_name]][[sample_size]][[type]][[prop_name]][[method]][[model]] <- NULL
            }
          }
        }
      }
    }
  }
}

# Combine into one tidy data frame for easier plotting/analysis
flatten_metrics_to_df <- function(metrics_list) {
  bind_rows(
    lapply(metrics_list, function(scen) {
      lapply(scen, function(sz) {
        lapply(sz, function(typ) {
          lapply(typ, function(prop) {
            lapply(prop, function(met) {
              lapply(met, function(mod) {
                mod # mod is a data frame or NULL
              }) %>% bind_rows()
            }) %>% bind_rows()
          }) %>% bind_rows()
        }) %>% bind_rows()
      }) %>% bind_rows()
    }) %>% bind_rows()
  )
}
metrics_df <- flatten_metrics_to_df(metrics)

saveRDS(metrics, paste0(res_path, "new_format/metrics_continuous_miss_nested.RDS"))
saveRDS(metrics_df, paste0(res_path, "new_format/metrics_continuous_miss_tidy.RDS"))
