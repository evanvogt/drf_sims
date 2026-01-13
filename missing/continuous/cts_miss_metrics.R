##########
# Title: evaluation metrics missing + cts outcome
##########

# libraries
library(here)
library(dplyr)

# paths
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live"

# parameters
scenarios <- paste0("scenario_", seq(1:5))
sample_sizes <- as.character(500)
types <- c("prognostic", "predictive", "both")
props <- as.character(0.3)
mechanisms <- c("MAR", "AUX")
methods <- c("complete_cases", "mean_imputation", "missforest", "regression", "missing_indicator", "IPW", "none")
models <- c("causal_forest", "dr_random_forest", "dr_oracle", "dr_semi_oracle", "dr_superlearner")

results <- readRDS(file.path(path, "results/new_format/missing_continuous_all.RDS"))

metrics <- data.frame(matrix(nrow = 0, ncol = 11))
colnames(metrics) <- c("scenario", "n", "prop", "type", "mechansim", "method", "model", "bias", "corr", "mse", "run")

for (scenario in scenarios) {

  for (n in sample_sizes) {

    for (type in types) {
      if (scenario == "scenario_1" & type != "prognostic") next
      
      for (prop in props) {

        for (mechanism in mechanisms) {

          for (method in methods) {
            sim_list <- results[[scenario]][[n]][[type]][[prop]][[mechanism]][[method]]
            
            runs <- names(sim_list)
            
            met_run <- data.frame(matrix(nrow = 0, ncol = 11))
            colnames(met_run) <- c("scenario", "n", "prop", "type", "mechansim", "method", "model", "bias", "corr", "mse", "run")
            
            for (run in runs) {
              sim_res <- sim_list[[run]]
              
              models_run <- intersect(names(sim_res), models)
              truth <- sim_res$truth
              
              met_models <- data.frame(matrix(nrow = 0, ncol = 11))
              colnames(met_models) <- c("scenario", "n", "prop", "type", "mechansim", "method", "model", "bias", "corr", "mse", "run")
              for (model in models_run) {
                model_tau <- sim_res[[model]]$tau
                true_tau <- truth$tau
                
                bias <- mean(true_tau - model_tau, na.rm = TRUE)
                corr <- ifelse(scenario != "scenario_1", cor(true_tau, model_tau, use = "pairwise.complete.obs"), 0)
                mse <- mean((true_tau - model_tau)^2, na.rm = TRUE)
                
                met <- data.frame(
                  scenario = scenario,
                  n = as.numeric(n),
                  prop = as.numeric(prop),
                  type = type,
                  mechanism = mechanism,
                  method = method,
                  model = model,
                  bias = bias,
                  corr = corr,
                  mse = mse,
                  run = run)
                met_models <- rbind(met_models, met)
              }
              met_run <- rbind(met_run, met_models)
            }
            metrics <- rbind(metrics, met_run)
          }
        }
      }
    }
  }
}

# tidy up the dataframe
metrics <- metrics %>%
  mutate(
    scenario = factor(scenario),
    type = factor(type),
    mechanism = factor(mechanism),
    method = factor(
      recode(method,
             complete_cases = "complete cases",
             mean_imputation = "mean imputation",
             missforest = "missForest",
             regression = "regression",
             missing_indicator = "missing indicator",
             IPW = "IPW",
             none = "in-built"),
      levels = c("complete cases", "mean imputation", "missForest", "regression", "missing indicator", "IPW", "in-built")
    ),
    model = factor(
      recode(model,
             causal_forest = "Causal forest",
             dr_random_forest = "DR-RandomForest",
             dr_superlearner = "DR-SuperLearner",
             dr_oracle = "DR-oracle",
             dr_semi_oracle = "DR-semi-oracle"),
      levels = c("Causal forest", "DR-RandomForest", "DR-SuperLearner", "DR-oracle", "DR-semi-oracle")
    ),
    run = as.numeric(run)
    )

# save metrics file
metric_outfile <- file.path(path, "results/new_format/metrics_cts_miss_df.RDS")
saveRDS(metrics, metric_outfile)
