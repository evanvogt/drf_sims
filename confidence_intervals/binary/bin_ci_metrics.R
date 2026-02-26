##########
# title: metrics for bin outcome confidence intervals
##########

# libraries
library(here)
library(dplyr)
library(tidyr)
library(purrr)

# paths
res_path <- file.path(dirname(here()), "results", "confidence_intervals", "binary")
out_file <- file.path(res_path, "ci_bin_metrics.RDS")

# parameters
models <- c("causal_forest", "dr_random_forest", "dr_oracle", "dr_semi_oracle")


# results list
all_results_df <- readRDS(file.path(res_path, "ci_bin_all.RDS"))

metrics <- all_results_df %>%
  # Unnest the simulation runs
  unnest_longer(results) %>%
  # one row per scenario / n / CI_sf / run
  mutate(run = map_int(results, ~.x$run),         # Extract run number
         sim_res = map(results, ~.x$result)) %>%  # Extract sim_res list
  select(-results) %>%
  # Map over models within each sim_res
  mutate(metrics = pmap(list(scenario, n, CI_sf, run, sim_res), function(scenario, n, CI_sf, run, sim_res) {
    models_run <- intersect(names(sim_res), models)
    truth <- sim_res$truth
    true_tau <- truth$tau
    # metrics by model
    map_dfr(models_run, function(model) {
      lower_bound <- sim_res[[model]]$hb_lb
      upper_bound <- sim_res[[model]]$hb_ub
      marginal_coverage <- mean(as.numeric(true_tau >= lower_bound & true_tau <= upper_bound))
      simultaneous_coverage <- as.numeric(all(true_tau >= lower_bound & true_tau <= upper_bound))
      mean_ci_length <- mean(upper_bound - lower_bound)
      met <- tibble(
        scenario = scenario,
        n = n,
        CI_sf = CI_sf,
        model = model,
        run = run,
        marginal_coverage = marginal_coverage,
        simultaneous_coverage = simultaneous_coverage,
        mean_ci_length = mean_ci_length
      )
      # Add causal forest inbuilt CIs
      if (model == "causal_forest" && !is.null(sim_res[[model]]$variance)) {
        var <- sim_res[[model]]$variance
        cf_tau <- sim_res[[model]]$tau
        lower_bound_cf <- cf_tau + qnorm(0.025) * sqrt(var)
        upper_bound_cf <- cf_tau + qnorm(0.975) * sqrt(var)
        marginal_coverage_cf <- mean(as.numeric(true_tau >= lower_bound_cf & true_tau <= upper_bound_cf))
        simultaneous_coverage_cf <- as.numeric(all(true_tau >= lower_bound_cf & true_tau <= upper_bound_cf))
        mean_ci_length <- mean(upper_bound_cf - lower_bound_cf)
        met <- bind_rows(
          met,
          tibble(
            scenario = scenario,
            n = n,
            CI_sf = CI_sf,
            model = "causal_forest_inbuilt",
            run = run,
            marginal_coverage = marginal_coverage_cf,
            simultaneous_coverage = simultaneous_coverage_cf,
            mean_ci_length = mean_ci_length
          )
        )
      }
      met
    })
  })) %>%
  # Unnest the metrics column to get final data frame
  select(metrics) %>%
  unnest(metrics)


# save metrics file
saveRDS(metrics, out_file)
print("metrics calculated!")