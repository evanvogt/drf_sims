##########
# title: metrics for single-event survival
##########
# As competing_risk/surv_metrics.R, on one target: every arm estimates the
# RMST CATE and is scored against tau_RMST, with cate_metrics() so the bias
# convention matches the rest of the repo.

library(here)
source(here("competing_risk", "single_event", "se_config.R"))
source(here("R", "metrics.R"))

# NOTE: frameworks_run below is intersect(names(sim_res), frameworks), so a
# framework missing from `frameworks` is dropped silently rather than erroring.
# Add a new arm here and to framework_info together.
framework_info <- tribble(
  ~framework,             ~learner, ~nuisance,
  "csf",                  "CSF",    "grf",
  "pseudo_cf_whole_oob",  "CF",     "RF",
  "pseudo_dr_whole_oob",  "DR",     "RF",
  "pseudo_t_whole_oob",   "T",      "RF",
  "sl_dr_whole",          "DR",     "SL",
  "sl_t_whole",           "T",      "SL",
  "rsf_dr_oob",           "DR",     "RSF",
  "rsf_t_oob",            "T",      "RSF",
  "rsf_dr_scf",           "DR",     "RSF (scf)",
  "rsf_t_scf",            "T",      "RSF (scf)"
)
frameworks <- framework_info$framework

all_results_df <- readRDS(file.path(study$res_path, "se_all.RDS"))

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

      true_tau       <- sim_res$truth$tau_RMST
      frameworks_run <- intersect(names(sim_res), frameworks)

      map_dfr(frameworks_run, function(framework) {
        model_tau <- sim_res[[framework]]$RMST
        bind_cols(
          tibble(scenario = scenario, n = n, censoring = censoring, run = run,
                 framework = framework, target = "RMST"),
          cate_metrics(model_tau, true_tau, scenario),
          # spread of the estimates: in scenario 1 (tau = 0) this is spurious
          # heterogeneity, the part of the RMSE that is not bias
          tibble(c_stat = c_statistic(model_tau, true_tau, scenario),
                 tau_sd = sd(model_tau, na.rm = TRUE))
        )
      })
    }
  )) %>%
  select(metrics) %>%
  unnest(metrics) %>%
  left_join(framework_info, by = "framework")

saveRDS(metrics, file.path(study$res_path, "se_metrics.RDS"))
print("metrics calculated!")
