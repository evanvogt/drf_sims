##########
# title: metrics for competing risk outcome
##########

# libraries
library(here)
library(dplyr)
library(tidyr)
library(purrr)

# paths
res_path <- file.path(dirname(here()), "results", "competing_risk")
out_file <- file.path(res_path, "surv_metrics.RDS")

# parameters
frameworks <- c("ipw", "csf_cs", "csf_sh", "pseudo_cf", "pseudo_dr")

# which targets are valid per framework
framework_targets <- list(
  ipw      = c("RMST1", "RMST2", "RMSTc"),
  csf_cs    = c("RMST1", "RMST2", "RMSTc"),
  csf_sh    = c("RMST1", "RMST2"),
  pseudo_cf = c("RMTL1", "RMTL2", "RMSTc"),
  pseudo_dr = c("RMTL1", "RMTL2", "RMSTc")
)

# framework-specific truth column mapping
# ipw and csf_cs remove competing events so they target the cause-specific (net) RMST = integral of S*(t)
# csf_sh keeps competing events in the risk set (Fine-Gray) so it targets the subdistribution RMST = horizon - RMTL
framework_truth_map <- list(
  ipw       = c(RMST1 = "tau_RMST1_cs", RMST2 = "tau_RMST2_cs", RMSTc = "tau_RMSTc"),
  csf_cs    = c(RMST1 = "tau_RMST1_cs", RMST2 = "tau_RMST2_cs", RMSTc = "tau_RMSTc"),
  csf_sh    = c(RMST1 = "tau_RMST1",    RMST2 = "tau_RMST2"),
  pseudo_cf = c(RMTL1 = "tau_RMTL1",    RMTL2 = "tau_RMTL2",    RMSTc = "tau_RMSTc"),
  pseudo_dr = c(RMTL1 = "tau_RMTL1",    RMTL2 = "tau_RMTL2",    RMSTc = "tau_RMSTc")
)

# results
all_results_df <- readRDS(file.path(res_path, "surv_all.RDS"))

metrics <- all_results_df %>%
  # unnest runs
  unnest_longer(results) %>%
  # one row per param-combo & run
  mutate(
    run     = map_int(results, ~.x$run),
    sim_res = map(results,     ~.x$result)
  ) %>%
  select(-results) %>%
  # map over frameworks x targets within each sim_res
  mutate(metrics = pmap(
    list(scenario, n, censoring, run, sim_res),
    function(scenario, n, censoring, run, sim_res) {
      
      truth          <- sim_res$truth
      frameworks_run <- intersect(names(sim_res), frameworks)
      
      map_dfr(frameworks_run, function(framework) {
        
        fw_data      <- sim_res[[framework]]
        targets_run  <- intersect(names(fw_data), framework_targets[[framework]])
        
        map_dfr(targets_run, function(target) {
          
          model_tau <- fw_data[[target]]
          true_tau  <- truth[[ framework_truth_map[[framework]][[target]] ]]
          
          bias <- mean(true_tau - model_tau, na.rm = TRUE)
          mse  <- mean((true_tau - model_tau)^2, na.rm = TRUE)
          corr <- ifelse(
            scenario != 1,
            cor(true_tau, model_tau, use = "pairwise.complete.obs"),
            0
          )
          
          tibble(
            scenario  = scenario,
            n         = n,
            censoring = censoring,
            run       = run,
            framework = framework,
            target    = target,
            bias      = bias,
            mse       = mse,
            corr      = corr
          )
        })
      })
    }
  )) %>%
  select(metrics) %>%
  unnest(metrics) %>%
  # combine RMST and RMTL metrics (inverses of eachother)
  mutate(target = case_when(target %in% c("RMST1", "RMTL1") ~ "Event 1",
                            target %in% c("RMST2", "RMTL2") ~ "Event 2",
                            target == "RMSTc" ~ "Combined"))

# save metrics file
saveRDS(metrics, out_file)
print("metrics calculated!")
