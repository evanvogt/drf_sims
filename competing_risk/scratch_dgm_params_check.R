
# Generate all scenarios and bind
scenarios <- 1:7
n <- 5000
censoring <- F
scenario_labels <- survival_scenario_params$description

data_all <- bind_rows(lapply(scenarios, function(s) {
  gen <- generate_surv_data(
    scenario = s,
    n = n,                 
    return_truth = FALSE,
    censoring = censoring   
  )
  df <- gen$dataset
  df$scenario <- factor(s, levels = scenarios, labels = scenario_labels)
  df
}))

# Tidy factors
plot_df <- data_all %>%
  mutate(
    D = factor(D, levels = c(0, 1, 2),
               labels = c("Censored", "Event 1", "Event 2")),
    W = factor(W, levels = c(0, 1),
               labels = c("Control", "Treated"))
  )

# Common bandwidth for comparability
bw_all <- stats::bw.nrd0(plot_df$Y)

# if no censoring, remove censored people
if (censoring == F) {
  plot_df <- plot_df %>%
    filter(D != "Censored") %>%
    droplevels()
}

# dist of event times
ggplot(plot_df, aes(x = Y, color = D, linetype = W)) +
  geom_density(bw = bw_all, n = 512, adjust = 1) +
  facet_wrap(~ scenario, ncol = 2, scales = "free_y") +
  labs(
    title = "Y density by scenario",
    x = "Time (Y)", y = "Density",
    color = "Event type", linetype = "Treatment"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# survival time
# ---------- tables of ATEs

scenarios <- 1:7
n <- 500
censoring <- FALSE

truth_tab <- bind_rows(lapply(scenarios, function(s) {
  gen <- generate_surv_data(
    scenario = s,
    n = n,
    return_truth = TRUE,
    censoring = censoring
  )
  tr <- gen$truth
  m <- nrow(tr)
  data.frame(
    scenario    = survival_scenario_params$description[survival_scenario_params$scenario == s],
    tau_RMSTc   = mean(tr$tau_RMSTc),
    se_RMSTc    = sd(tr$tau_RMSTc) / sqrt(m),
    lcl_RMSTc   = mean(tr$tau_RMSTc) - 1.96 * sd(tr$tau_RMSTc) / sqrt(m),
    ucl_RMSTc   = mean(tr$tau_RMSTc) + 1.96 * sd(tr$tau_RMSTc) / sqrt(m),
    
    tau_RMST1   = mean(tr$tau_RMST1),
    se_RMST1    = sd(tr$tau_RMST1) / sqrt(m),
    lcl_RMST1   = mean(tr$tau_RMST1) - 1.96 * sd(tr$tau_RMST1) / sqrt(m),
    ucl_RMST1   = mean(tr$tau_RMST1) + 1.96 * sd(tr$tau_RMST1) / sqrt(m),
    
    tau_RMST2   = mean(tr$tau_RMST2),
    se_RMST2    = sd(tr$tau_RMST2) / sqrt(m),
    lcl_RMST2   = mean(tr$tau_RMST2) - 1.96 * sd(tr$tau_RMST2) / sqrt(m),
    ucl_RMST2   = mean(tr$tau_RMST2) + 1.96 * sd(tr$tau_RMST2) / sqrt(m),
    
    tau_sh      = mean(tr$tau_sh),
    se_sh       = sd(tr$tau_sh) / sqrt(m),
    lcl_sh      = mean(tr$tau_sh) - 1.96 * sd(tr$tau_sh) / sqrt(m),
    ucl_sh      = mean(tr$tau_sh) + 1.96 * sd(tr$tau_sh) / sqrt(m)
  )
}))

numeric_cols <- setdiff(colnames(truth_tab), "scenario")
truth_tab <- truth_tab %>%
  mutate_at(numeric_cols, function(x) {signif(x, 2)})

truth_tab %>% select(contains(c("scenario", "tau"))) %>% View()
truth_tab
