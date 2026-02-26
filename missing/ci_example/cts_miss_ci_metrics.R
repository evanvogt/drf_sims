##########
# Title: evaluation metrics missing + cts outcome
##########

# libraries
library(here)
library(dplyr)
library(tidyr)
library(purrr)

# paths
res_path <- file.path(dirname(here()), "results", "missing", "ci_example")
out_file <- file.path(res_path, "cts_miss_ci_metrics.RDS")

# parameters
models <- c("causal_forest", "dr_random_forest", "dr_oracle", "dr_semi_oracle", "dr_superlearner")

# results list
all_results_df <- readRDS(file.path(res_path, "cts_miss_ci_all.RDS"))

metrics <- all_results_df %>%
  # unnest runs
  unnest_longer(results) %>%
  # one row per param-combo & run
  mutate(run = map_int(results, ~.x$run),         # Extract run number
         sim_res = map(results, ~.x$result)) %>%  # Extract sim_res list
  select(-results) %>%
  # Map over models within each sim_res
  mutate(metrics = pmap(list(scenario, n, type, prop, mechanism, method, run, sim_res), function(scenario, n, type, prop, mechanism, method, run, sim_res) {
    models_run <- intersect(names(sim_res), models)
    truth <- sim_res$truth
    true_tau <- truth$tau
    # metrics by model
    map_dfr(models_run, function(model) {
      # CATE performance
      model_tau <- sim_res[[model]]$tau
      bias <- mean(true_tau - model_tau, na.rm = T)
      mse <- mean((true_tau - model_tau)^2, na.rm = T)

      # Confidence Intervals
      # Pooled sample method
      lb_p <- sim_res[[model]]$lb_pooled
      ub_p <- sim_res[[model]]$ub_pooled
      mar_cov_p <- mean(as.numeric(true_tau >= lb_p & true_tau <= ub_p))
      sim_cov_p <- as.numeric(all(true_tau >= lb_p & true_tau <= ub_p))
      m_ci_len_p <- mean(ub_p - lb_p)
      
      # MI + boot method
      lb_mib <- sim_res[[model]]$lb_mib
      ub_mib <- sim_res[[model]]$ub_mib
      mar_cov_mib <- mean(as.numeric(true_tau >= lb_mib & true_tau <= ub_mib))
      sim_cov_mib <- as.numeric(all(true_tau >= lb_mib & true_tau <= ub_mib))
      m_ci_len_mib <- mean(ub_mib - lb_mib)
      
      # Hybrid method
      lb_h <- sim_res[[model]]$lb_hybrid
      ub_h <- sim_res[[model]]$ub_hybrid
      mar_cov_h <- mean(as.numeric(true_tau >= lb_h & true_tau <= ub_h))
      sim_cov_h <- as.numeric(all(true_tau >= lb_h & true_tau <= ub_h))
      m_ci_len_h <- mean(ub_h - lb_h)
      
      tibble(
        scenario = scenario,
        n = n,
        type = type,
        prop = prop,
        mechanism = mechanism,
        method = method,
        model = model,
        run = run,
        bias = bias,
        mse = mse,
        CI_method = c("pooled", "MI boot", "hybrid"),
        marginal_coverage = c(mar_cov_p, mar_cov_mib, mar_cov_h),
        simultaneous_coverage = c(sim_cov_p, sim_cov_mib, sim_cov_h),
        mean_ci_length = c(m_ci_len_p, m_ci_len_mib, m_ci_len_h)
      )
    })
  })) %>% select(metrics) %>%
  unnest(metrics)

# save metrics file
saveRDS(metrics, out_file)
print("metrics calculated!")
