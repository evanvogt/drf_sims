##########
# title: figures for the thesis chapter - cts outcome confidence intervals
##########

# libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(paletteer)
library(here)
library(patchwork)
library(purrr)
library(ggridges)
library(scales)

# paths
path <- here()
res_path <- file.path(dirname(path), "results", "confidence_intervals", "continuous")
fig_path <- file.path(dirname(path), "results", "thesis_figures", "cts_ci")
dir.create(fig_path, showWarnings = F, recursive = T)

metrics <- readRDS(file.path(res_path, "ci_cts_metrics.RDS"))

metrics <- metrics %>%
  filter(scenario %in% c(1, 3, 8, 9)) %>%
  mutate(
    scenario = factor(
      case_when(scenario == 1 ~ "Null",
                scenario == 3 ~ "Simple",
                scenario == 8 ~ "Complex",
                scenario == 9 ~ "Non-linear"), levels = c("Null", "Simple", "Complex", "Non-linear")),
    n = factor(n, levels = c(500, 1000)),
    CI_sf = factor(CI_sf, levels = seq(0.05, 0.5, 0.05)),
    model = factor(
      recode(model,
             causal_forest = "Causal forest",
             causal_forest_inbuilt = "Causal forest (inbuilt CIs)",
             dr_random_forest = "DR-RandomForest",
             dr_oracle = "DR-oracle",
             dr_semi_oracle = "DR-semi-oracle"),
      levels = c("Causal forest", "Causal forest (inbuilt CIs)", "DR-RandomForest", "DR-oracle", "DR-semi-oracle")
    )
  )

metrics_summary <- metrics %>%
  group_by(scenario, n, CI_sf, model) %>%
  summarise(
    mean_marg_cov = mean(marginal_coverage, na.rm = T),
    mcse_marg_cov = sd(marginal_coverage, na.rm = T)/sqrt(n()),
    mean_simul_cov = mean(simultaneous_coverage, na.rm = T),
    mcse_simul_cov = sd(simultaneous_coverage, na.rm = T)/sqrt(n()),
    mean_ci_len = mean(mean_ci_length, na.rm = T),
    mcse_ci_len = sd(mean_ci_length, na.rm = T)/sqrt(n()),
    .groups = "drop"
  )


mc_sum_plot <- metrics_summary %>%
  ggplot(aes(x = CI_sf, y = mean_marg_cov, colour = model, ymin = mean_marg_cov - mcse_marg_cov, ymax = mean_marg_cov + mcse_marg_cov)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(position = position_dodge(width = 0.5)) +
  facet_grid(cols = vars(scenario), rows = vars(n), scales = "free_y", axes = "all_x") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "Marginal coverage",
       x = "Subsampling ratio") +
  theme(axis.text.x = element_text(size = 4.5))
ggsave("cts_ci_mar_cov_sum.png", path = fig_path, width = 21, height = 15, units = "cm")

sc_sum_plot <- metrics_summary %>%
  ggplot(aes(x = CI_sf, y = mean_simul_cov, colour = model, ymin = mean_simul_cov - mcse_simul_cov, ymax = mean_simul_cov + mcse_simul_cov)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(position = position_dodge(width = 0.5)) +
  facet_grid(cols = vars(scenario), rows = vars(n), scales = "free_y", axes = "all_x") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "Simultaneous coverage",
       x = "Subsampling ratio")
ggsave("cts_ci_simul_cov_sum.png", path = fig_path, width = 21, height = 15, units = "cm")

cl_sum_plot <- metrics_summary %>%
  ggplot(aes(x = CI_sf, y = mean_ci_len, colour = model, ymin = mean_ci_len - mcse_ci_len, ymax = mean_ci_len + mcse_ci_len)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(position = position_dodge(width = 0.5)) +
  facet_grid(cols = vars(scenario), rows = vars(n), scales = "free_y", axes = "all_x") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "Confidence interval length",
       x = "Subsampling ratio")
ggsave("cts_ci_ci_len_sum.png", path = fig_path, width = 21, height = 15, units = "cm")
