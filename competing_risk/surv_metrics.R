##########
# title: metrics for competing risk outcome
##########
# This study does not use compute_metrics(): its results nest as
# framework x target rather than one entry per model, and each target is scored
# against a different truth column. The point metrics themselves come from
# cate_metrics() so the bias convention matches the rest of the repo.

library(here)
source(here("competing_risk/surv_config.R"))
source(here("R", "metrics.R"))

frameworks <- c("ipw", "csf_cs", "csf_sh", "pseudo_cf", "pseudo_dr")

# which targets are valid per framework
framework_targets <- list(
  ipw       = c("RMST1", "RMST2", "RMSTc"),
  csf_cs    = c("RMST1", "RMST2", "RMSTc"),
  csf_sh    = c("RMST1", "RMST2"),
  pseudo_cf = c("RMTL1", "RMTL2", "RMSTc"),
  pseudo_dr = c("RMTL1", "RMTL2", "RMSTc")
)

# framework-specific truth column mapping.
# ipw and csf_cs remove competing events, so they target the cause-specific (net)
# RMST = integral of S*(t). csf_sh keeps competing events in the risk set
# (Fine-Gray) so it targets the subdistribution RMST = horizon - RMTL.
framework_truth_map <- list(
  ipw       = c(RMST1 = "tau_RMST1_cs", RMST2 = "tau_RMST2_cs", RMSTc = "tau_RMSTc"),
  csf_cs    = c(RMST1 = "tau_RMST1_cs", RMST2 = "tau_RMST2_cs", RMSTc = "tau_RMSTc"),
  csf_sh    = c(RMST1 = "tau_RMST1",    RMST2 = "tau_RMST2"),
  pseudo_cf = c(RMTL1 = "tau_RMTL1",    RMTL2 = "tau_RMTL2",    RMSTc = "tau_RMSTc"),
  pseudo_dr = c(RMTL1 = "tau_RMTL1",    RMTL2 = "tau_RMTL2",    RMSTc = "tau_RMSTc")
)

all_results_df <- readRDS(file.path(study$res_path, "surv_all.RDS"))

# C-statistic = (Kendall tau_b + 1) / 2, equivalent to Harrell's C for a
# continuous outcome. Undefined in the null scenario, where it is 0.5.
c_statistic <- function(est, true, scenario) {
  if (scenario == 1) return(0.5)
  (cor(true, est, method = "kendall", use = "pairwise.complete.obs") + 1) / 2
}

metrics <- all_results_df %>%
  unnest_longer(results) %>%
  mutate(run     = map_int(results, ~ .x$run),
         sim_res = map(results,     ~ .x$result)) %>%
  select(-results) %>%
  mutate(metrics = pmap(
    list(scenario, n, censoring, run, sim_res),
    function(scenario, n, censoring, run, sim_res) {

      truth          <- sim_res$truth
      frameworks_run <- intersect(names(sim_res), frameworks)

      map_dfr(frameworks_run, function(framework) {

        fw_data     <- sim_res[[framework]]
        targets_run <- intersect(names(fw_data), framework_targets[[framework]])

        map_dfr(targets_run, function(target) {

          model_tau <- fw_data[[target]]
          true_tau  <- truth[[framework_truth_map[[framework]][[target]]]]

          bind_cols(
            tibble(scenario = scenario, n = n, censoring = censoring, run = run,
                   framework = framework, target = target),
            cate_metrics(model_tau, true_tau, scenario),
            tibble(c_stat = c_statistic(model_tau, true_tau, scenario))
          )
        })
      })
    }
  )) %>%
  select(metrics) %>%
  unnest(metrics) %>%
  # RMST and RMTL are inverses of each other, so label by event rather than scale
  mutate(target = case_when(target %in% c("RMST1", "RMTL1") ~ "Event 1",
                            target %in% c("RMST2", "RMTL2") ~ "Event 2",
                            target == "RMSTc" ~ "Combined"))

saveRDS(metrics, file.path(study$res_path, "surv_metrics.RDS"))
print("metrics calculated!")
