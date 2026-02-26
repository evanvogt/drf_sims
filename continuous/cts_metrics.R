##########
# title: metrics for cts outcome
##########

# libraries
library(here)
library(dplyr)
library(tidyr)
library(purrr)

# paths
res_path <- file.path(dirname(here()), "results", "continuous")
out_file <- file.path(res_path, "cts_metrics.RDS")

# parameters
models <- c("causal_forest", "dr_random_forest", "dr_oracle", "dr_semi_oracle", "dr_superlearner")

# results list
all_results_df <- readRDS(file.path(res_path, "cts_all.RDS"))

metrics <- all_results_df %>%
  # unnest runs
  unnest_longer(results) %>%
  # one row per param-combo & run
  mutate(run = map_int(results, ~.x$run),         # Extract run number
         sim_res = map(results, ~.x$result)) %>%  # Extract sim_res list
  select(-results) %>%
  # map over models within each sim_res
  mutate(metrics = pmap(list(scenario, n, run, sim_res), function(scenario, n, run, sim_res) {
    models_run <- intersect(names(sim_res), models)
    truth <- sim_res$truth
    true_tau <- truth$tau
    # metrics by model
    map_dfr(models_run, function(model) {
      # CATE performance
      model_tau <- sim_res[[model]]$tau
      bias <- mean(true_tau - model_tau, na.rm = T)
      mse <- mean((true_tau - model_tau)^2, na.rm = T)
      corr <- ifelse(scenario != 1, cor(true_tau, model_tau, use = "pairwise.complete.obs"), 0)
      
      # HTE test metrics - add in once sims have rerun with HTE tests
      BLP_p <- sim_res[[model]]$BLP_whole[4,2]
      indep_cate <- sim_res[[model]]$independence_cate$p_value
      indep_po <- sim_res[[model]]$independence_po$p_value
      tibble(
        scenario = scenario,
        n = n,
        model = model,
        run = run,
        bias = bias,
        mse = mse,
        corr = corr,
        BLP_p = BLP_p,
        indep_cate = indep_cate,
        indep_po = indep_po
      )
    })
  })) %>% select(metrics) %>%
  unnest(metrics)

# save metrics file
saveRDS(metrics, out_file)
print("metrics calculated!")