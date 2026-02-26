##########
# title: figures for the thesis chapter - cts miss CI
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
res_path <- file.path(dirname(path), "results", "missing", "ci_example")
fig_path <- file.path(dirname(path), "results", "thesis_figures", "cts_ss")


metrics <- readRDS(file.path(res_path, "cts_miss_ci_metrics.RDS"))

# tidy up
metrics <- metrics %>%
  filter(scenario %in% c(1, 2, 4, 5)) %>%
  mutate(scenario = factor(
    case_when(scenario == 1 ~ "Null",
              scenario == 2 ~ "Simple",
              scenario == 4 ~ "Complex",
              scenario == 5 ~ "Non-linear"), levels = c("Null", "Simple", "Complex", "Non-linear")),
    n = factor(n, levels = c(100, 250, 500, 1000)),
    model = factor(
      recode(model,
             causal_forest = "Causal forest",
             dr_random_forest = "DR-RandomForest",
             dr_oracle = "DR-oracle",
             dr_semi_oracle = "DR-semi-oracle"),
      levels = c("Causal forest", "DR-RandomForest", "DR-oracle", "DR-semi-oracle")),
    CI_method = factor(CI_method, levels = c("pooled", "MI boot", "hybrid")))


# per scenario summaries
metrics_summary <- metrics %>%
  group_by(scenario, n, type, prop, mechanism, method, model, CI_method) %>%
  summarise(
    mean_bias = mean(bias, na.rm = T),
    mcse_bias = sd(bias, na.rm = T)/sqrt(n()),
    mean_mse = mean(mse, na.rm = T),
    mcse_mse = sd(mse, na.rm = T)/sqrt(n()),
    mean_mar_cov = mean(marginal_coverage, na.rm = T),
    mcse_mar_cov = sd(marginal_coverage, na.rm = T),
    mean_sim_cov = mean(simultaneous_coverage, na.rm = T),
    mcse_sim_cov = sd(simultaneous_coverage, na.rm = T),
    mean_ci_len = mean(mean_ci_length, na.rm = T),
    mcse_ci_len = sd(mean_ci_length, na.rm = T),
    .groups = "drop"
  )

# combining everything into a single plot?


# marginal coverage
mar_cov_plot <- metrics %>%
  ggplot(aes(x=CI_method, y = marginal_coverage, colour = model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_wrap(~scenario) +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "Marginal coverage",
       x = "Method")

# summary
mar_cov_sum_plot <- metrics_summary %>%
  ggplot(aes(x = CI_method, y = mean_mar_cov, colour = model, ymin = mean_mar_cov - mcse_mar_cov, ymax = mean_mar_cov + mcse_mar_cov)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(position = position_dodge(width = 0.5), linewidth = 0.3) +
  facet_wrap(~scenario) +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "Marginal Coverage",
       x = "Method")

# simultaneous coverage
sim_cov_plot <- metrics %>%
  ggplot(aes(x=CI_method, y = simultaneous_coverage, colour = model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_wrap(~scenario) +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "Simultaneous coverage",
       x = "Method") +
  theme(axis.text.x = element_text(angle = 90))

# summary
sim_cov_sum_plot <- metrics_summary %>%
  ggplot(aes(x = CI_method, y = mean_sim_cov, colour = model, ymin = mean_sim_cov - mcse_sim_cov, ymax = mean_sim_cov + mcse_sim_cov)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(position = position_dodge(width = 0.5), linewidth = 0.3) +
  facet_wrap(~scenario) +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "Simultaneous coverage",
       x = "Method")

# CI length
ci_len_plot <- metrics %>%
  ggplot(aes(x=CI_method, y = mean_ci_length, colour = model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_wrap(~scenario) +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "CI length",
       x = "Method")

# summary
ci_len_sum_plot <- metrics_summary %>%
  ggplot(aes(x = CI_method, y = mean_ci_len, colour = model, ymin = mean_ci_len - mcse_ci_len, ymax = mean_ci_len + mcse_ci_len)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(position = position_dodge(width = 0.5), linewidth = 0.3) +
  facet_wrap(~scenario) +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "CI length",
       x = "Method")

# combine the plots
mar <- mar_cov_sum_plot + facet_grid(cols = vars(scenario)) + theme(axis.text.x = element_text(angle = 90)) + scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8))
sim <- sim_cov_sum_plot + facet_grid(cols = vars(scenario)) + theme(axis.text.x = element_text(angle = 90)) + ylim(0,1)
len <- ci_len_sum_plot + facet_grid(cols = vars(scenario)) + theme(axis.text.x = element_text(angle = 90)) 

(mar / len) + plot_layout(guides = "collect")
ggsave("miss_cts_ci_all.png", path = fig_path, width = 21, height = 15, units = "cm")
