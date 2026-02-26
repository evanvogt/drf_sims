##########
# title: figures for the thesis chapter - cts outcomes
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
res_path <- file.path(dirname(path), "results", "continuous")
fig_path <- file.path(dirname(path), "results", "thesis_figures", "cts_ss")


metrics <- readRDS(file.path(res_path, "cts_metrics.RDS"))

# tidy up
metrics <- metrics %>%
  filter(scenario %in% c(1, 3, 8, 9)) %>%
  mutate(scenario = factor(
    case_when(scenario == 1 ~ "Null",
              scenario == 3 ~ "Simple",
              scenario == 8 ~ "Complex",
              scenario == 9 ~ "Non-linear"), levels = c("Null", "Simple", "Complex", "Non-linear")),
    n = factor(n, levels = c(100, 250, 500, 1000)),
    model = factor(
      recode(model,
             causal_forest = "Causal forest",
             dr_random_forest = "DR-RandomForest",
             dr_oracle = "DR-oracle",
             dr_semi_oracle = "DR-semi-oracle",
             dr_superlearner = "DR-SuperLearner"),
      levels = c("Causal forest", "DR-RandomForest", "DR-oracle", "DR-semi-oracle", "DR-SuperLearner")))


# per scenario summaries
metrics_summary <- metrics %>%
  group_by(scenario, n, model) %>%
  summarise(
    mean_bias = mean(bias, na.rm = T),
    mcse_bias = sd(bias, na.rm = T)/sqrt(n()),
    mean_mse = mean(mse, na.rm = T),
    mcse_mse = sd(mse, na.rm = T)/sqrt(n()),
    mean_corr = mean(corr, na.rm = T),
    mcse_corr = sd(corr, na.rm = T)/sqrt(n()),
    mean_BLP = mean(BLP_p, na.rm = T),
    mcse_BLP = sd(BLP_p, na.rm = T)/sqrt(n()),
    mean_indep_cate = mean(indep_cate, na.rm = T),
    mcse_indep_cate = sd(indep_cate, na.rm = T)/sqrt(n()),
    mean_indep_po = mean(indep_po, na.rm = T),
    mcse_indep_po = sd(indep_po, na.rm = T)/sqrt(n()),
    .groups = "drop"
  )

# bias plots
bias_plot <- metrics %>%
  ggplot(aes(x=n, y = bias, colour = model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_wrap(~scenario) +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(title = "Average Bias distribution in CATE estimates",
       y = "Average bias",
       x = "Sample size")
ggsave("cts_bias_all.png", path = fig_path, width = 21, height = 15, units = "cm")

bias_sum_plot <- metrics_summary %>%
  ggplot(aes(x = n, y = mean_bias, colour = model, ymin = mean_bias - mcse_bias, ymax = mean_bias + mcse_bias)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = mean_bias - mcse_bias,
                    ymax = mean_bias + mcse_bias),
                position = position_dodge(width = 0.5), linewidth = 0.3) +
  facet_wrap(~scenario) +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(title = "Average Bias distribution in CATE estimates",
       y = "Average bias",
       x = "Sample size")
ggsave("cts_bias_summary.png", path = fig_path, width = 21, height = 15, units = "cm")

# MSE plots
# not super useful because the superlearner is so bad at n= 100
mse_plot <- metrics %>%
  filter(mse < 20) %>%
  ggplot(aes(x=n, y = mse, colour = model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_wrap(~scenario) +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(title = "Average MSE distribution in CATE estimates",
       y = "Average MSE",
       x = "Sample size")
ggsave("cts_mse_all.png", path = fig_path, width = 21, height = 15, units = "cm")

mse_sum_plot <- metrics_summary %>%
  #filter(!(n == 100 & model == "DR-SuperLearner")) %>%
  ggplot(aes(x = n, y = mean_mse, colour = model, ymin = mean_mse - mcse_mse, ymax = mean_mse + mcse_mse)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(position = position_dodge(width = 0.5), linewidth = 0.3) +
  facet_wrap(~scenario) +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(title = "Average MSE distribution in CATE estimates",
       y = "Average MSE",
       x = "Sample size")
ggsave("cts_mse_summary.png", path = fig_path, width = 21, height = 15, units = "cm")

#combining bias and mse into one plot
new_bias_plot <- bias_sum_plot + facet_grid(cols = vars(scenario)) + labs(title = "")
new_mse_plot <- mse_sum_plot + facet_grid(cols = vars(scenario)) + labs(title = "")

(new_bias_plot / new_mse_plot) + plot_layout(guides = "collect")
ggsave("cts_bias_mse.png", path = fig_path, width = 21, height = 15, units = "cm")

# correlation
corr_plot <- metrics %>%
  filter(scenario != "Null") %>%
  ggplot(aes(x = n, y = corr, colour = model)) +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  facet_wrap(~scenario) +
  theme_minimal() +
  labs(title = "Correlations of true vs estimated CATEs",
       y = "Correlation",
       x = "Sample size")
ggsave("cts_corr_all.png", path = fig_path, width = 21, height = 15, units = "cm")


corr_sum_plot <- metrics_summary %>%
  filter(scenario != "Null") %>%
  ggplot(aes(x = n, y = mean_corr, colour = model)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = mean_corr - mcse_corr,
                    ymax = mean_corr + mcse_corr),
                position = position_dodge(width = 0.5), linewidth = 0.3) +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  facet_wrap(~scenario) +
  theme_minimal() +
  labs(title = "Mean correlations of true vs estiated CATEs",
       y = "Correlation",
       x = "sample_size")
ggsave("cts_corr_summary.png", path = fig_path, width = 21, height = 15, units = "cm")

# BLP test
BLP_plot <- metrics %>%
  ggplot(aes(x = BLP_p, y = model, color = model)) +
  geom_density_ridges(alpha = 0.5, fill = "transparent") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  facet_grid(rows = vars(n), cols = vars(scenario)) +
  theme_minimal() +
  labs(title = "Distributions of BLP test p-values",
       y = "Model",
       x = "BLP p-value") +
  theme(legend.position = "none")
ggsave("cts_blp_all.png", path = fig_path, width = 21, height = 15, units = "cm")

BLP_sum_plot <- metrics_summary %>%
  ggplot(aes(x = n, y = mean_BLP, colour = model)) +
  geom_hline(yintercept = 0.1, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = mean_BLP - mcse_BLP,
                    ymax = mean_BLP + mcse_BLP),
                position = position_dodge(width = 0.5), linewidth = 0.3) +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  facet_wrap(~scenario) +
  theme_minimal() +
  labs(title = "BLP test p-values",
       y = "p-value",
       x = "sample_size")
ggsave("cts_blp_summary.png", path = fig_path, width = 21, height = 15, units = "cm")

# independence test on the CATE
indep_cate_plot <- metrics %>%
  ggplot(aes(x = indep_cate, y = model, color = model)) +
  geom_density_ridges(alpha = 0.5, fill = "transparent") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  facet_grid(rows = vars(n), cols = vars(scenario)) +
  theme_minimal() +
  labs(title = "Distributions of CATE-independence test p-values",
       y = "Model",
       x = "Independence test p-value") +
  theme(legend.position = "none") +
  xlim(0,1) +
ggsave("cts_indep_cate_all.png", path = fig_path, width = 21, height = 15, units = "cm")

indep_cate_sum_plot <- metrics_summary %>%
  ggplot(aes(x = n, y = mean_indep_cate, colour = model)) +
  geom_hline(yintercept = 0.1, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = mean_indep_cate - mcse_indep_cate,
                    ymax = mean_indep_cate + mcse_indep_cate),
                position = position_dodge(width = 0.5), linewidth = 0.3) +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  facet_wrap(~scenario) +
  theme_minimal() +
  labs(title = "CATE-independence test p-values",
       y = "p-value",
       x = "sample_size")
ggsave("cts_indep_cate_summary.png", path = fig_path, width = 21, height = 15, units = "cm")

# independence test on the potential outcomes
indep_po_plot <- metrics %>%
  ggplot(aes(x = indep_po, y = model, color = model)) +
  geom_density_ridges(alpha = 0.5, fill = "transparent") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  facet_grid(rows = vars(n), cols = vars(scenario)) +
  theme_minimal() +
  labs(title = "Distributions of PO-independence test p-values",
       y = "Model",
       x = "Independence test p-value") +
  xlim(0,1) +
  theme(legend.position = "none")
ggsave("cts_indep_po_all.png", path = fig_path, width = 21, height = 15, units = "cm")

indep_po_sum_plot <- metrics_summary %>%
  ggplot(aes(x = n, y = mean_indep_po, colour = model)) +
  geom_hline(yintercept = 0.1, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = mean_indep_po - mcse_indep_po,
                    ymax = mean_indep_po + mcse_indep_po),
                position = position_dodge(width = 0.5), linewidth = 0.3) +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  facet_wrap(~scenario) +
  theme_minimal() +
  labs(title = "PO-independence test p-values",
       y = "p-value",
       x = "sample_size")
ggsave("cts_indep_po_summary.png", path = fig_path, width = 21, height = 15, units = "cm")

# trying out plotting all of the different types of test on one dataset
all_tests_plot <- metrics %>%
  mutate(indep_cate = as.numeric(indep_cate),
         indep_po = as.numeric(indep_po)) %>%
  pivot_longer(cols = c(BLP_p, indep_cate, indep_po),
               names_to = "test",
               values_to = "pvalue") %>%
  ggplot(aes(x = pvalue, y = test, colour = model)) +
  geom_density_ridges(alpha = 0.1, fill = "transparent") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  facet_grid(rows = vars(n), cols = vars(scenario)) +
  theme_minimal() +
  labs(y = "Model",
       x = "Independence test p-value") +
  xlim(0,1)
ggsave("cts_tests_all.png", path = fig_path, width = 21, height = 15, units = "cm")

all_tests_sum_plot <- metrics_summary %>%
  pivot_longer(
    cols = c(mean_BLP, mcse_BLP, mean_indep_cate, mcse_indep_cate, mean_indep_po, mcse_indep_po),
    names_to = c("stat_type", "variable"),
    names_pattern = "(mean|mcse)_(.+)",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = stat_type,
    values_from = value
  ) %>%
  ggplot(aes(x = model, y = mean, colour = variable, ymin = mean - mcse, ymax = mean + mcse)) +
  geom_hline(yintercept = 0.1, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(position = position_dodge(width = 0.5), linewidth = 0.3) +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  facet_grid(rows = vars(n), cols = vars(scenario)) +
  theme_classic() +
  labs(title = "PO-independence test p-values",
       y = "p-value",
       x = "sample_size") +
  theme(axis.text.x = element_text(angle = 90))

  
