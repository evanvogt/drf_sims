##########
# title: figures for the continuous outcome study
##########
rm(list = ls())
# libraries
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(paletteer)
library(here)
library(patchwork)
library(purrr)

# paths
path <- here()
res_path <- file.path(dirname(path), "results", "continuous")
fig_path <- file.path(dirname(path), "results", "thesis_figures", "continuous")
dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)

metrics <- readRDS(file.path(res_path, "cts_metrics.RDS"))

# display order and labels - continuous has no shared label-dictionary script
# (unlike crossfitting/cf_metrics.R, which does double duty as pipeline +
# labels), so these stay inline, matching cts_results.Rmd's precedent
scenario_labels <- c(
  `1` = "1: No HTE",
  `2` = "2: Binary HTE",
  `3` = "3: Continuous HTE",
  `4` = "4: Two HTE vars",
  `5` = "5: Cts x binary interaction",
  `6` = "6: Two vars + cts x binary",
  `7` = "7: Cts x cts interaction",
  `8` = "8: Two vars + cts x cts",
  `9` = "9: Cosine HTE",
  `10` = "10: Exponential HTE"
)

model_labels <- c(
  causal_forest = "Causal forest",
  dr_random_forest = "DR-RandomForest",
  dr_oracle = "DR-oracle",
  dr_semi_oracle = "DR-semi-oracle",
  dr_superlearner = "DR-SuperLearner"
)

# tidy up
metrics <- metrics %>%
  mutate(
    scenario = factor(scenario, levels = names(scenario_labels) %>% as.integer(),
                      labels = unname(scenario_labels)),
    n = factor(n, levels = c(100, 250, 500, 1000)),
    model = factor(recode(model, !!!model_labels), levels = unname(model_labels))
  ) %>%
  droplevels()

# per (scenario, model, n) summaries
metrics_summary <- metrics %>%
  group_by(scenario, model, n) %>%
  summarise(
    mean_bias = mean(bias, na.rm = T),
    mcse_bias = sd(bias, na.rm = T) / sqrt(n()),
    mean_mse = mean(mse, na.rm = T),
    mcse_mse = sd(mse, na.rm = T) / sqrt(n()),
    mean_corr = mean(corr, na.rm = T),
    mcse_corr = sd(corr, na.rm = T) / sqrt(n()),
    mean_BLP = mean(BLP_p, na.rm = T),
    mcse_BLP = sd(BLP_p, na.rm = T) / sqrt(n()),
    mean_indep_cate = mean(indep_cate, na.rm = T),
    mcse_indep_cate = sd(indep_cate, na.rm = T) / sqrt(n()),
    mean_indep_po = mean(indep_po, na.rm = T),
    mcse_indep_po = sd(indep_po, na.rm = T) / sqrt(n()),
    .groups = "drop"
  )

# helper: the house summary plot, points with MCSE error bars, faceted by
# scenario (there's no family/variant dimension here - just scenario x model
# x n) - see crossfitting/cf_results.R's summary_plot for the pattern this
# follows
summary_plot <- function(df, mean_col, mcse_col, title, ylab, hline = 0) {
  df %>%
    mutate(est = .data[[mean_col]],
           lo = .data[[mean_col]] - .data[[mcse_col]],
           hi = .data[[mean_col]] + .data[[mcse_col]]) %>%
    ggplot(aes(x = n, y = est, colour = model, ymin = lo, ymax = hi)) +
    geom_hline(yintercept = hline, linetype = "dashed") +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbar(position = position_dodge(width = 0.5), linewidth = 0.3, width = 0.3) +
    facet_wrap(~scenario, scales = "free") +
    scale_colour_paletteer_d("rcartocolor::Safe") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title, y = ylab, x = "Sample size", colour = "Model")
}

# --- bias --------------------------------------------------------------
bias_plot <- metrics %>%
  ggplot(aes(x = n, y = bias, colour = model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_wrap(~scenario, scales = "free") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Bias in CATE estimates", y = "Bias", x = "Sample size", colour = "Model")
ggsave("cts_bias_all.png", path = fig_path, width = 21, height = 15, units = "cm")

bias_sum_plot <- summary_plot(metrics_summary, "mean_bias", "mcse_bias",
                              "Mean bias in CATE estimates", "Mean bias")
ggsave("cts_bias_summary.png", plot = bias_sum_plot, path = fig_path,
       width = 21, height = 15, units = "cm")

# --- MSE -----------------------------------------------------------------
mse_plot <- metrics %>%
  ggplot(aes(x = n, y = mse, colour = model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_wrap(~scenario, scales = "free") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "MSE of CATE estimates", y = "MSE", x = "Sample size", colour = "Model")
ggsave("cts_mse_all.png", path = fig_path, width = 21, height = 15, units = "cm")

mse_sum_plot <- summary_plot(metrics_summary, "mean_mse", "mcse_mse",
                             "Mean MSE of CATE estimates", "Mean MSE")
ggsave("cts_mse_summary.png", plot = mse_sum_plot, path = fig_path,
       width = 21, height = 15, units = "cm")

# --- correlation with the truth ---------------------------------------------
# scenario 1 has no heterogeneity, so corr is undefined and forced to 0 -
# excluded here, same convention as crossfitting's results files
corr_df <- filter(metrics, scenario != scenario_labels[["1"]])

corr_plot <- corr_df %>%
  ggplot(aes(x = n, y = corr, colour = model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_wrap(~scenario, scales = "free_y") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation between estimated and true CATEs",
       y = "Pearson correlation", x = "Sample size", colour = "Model")
ggsave("cts_corr_all.png", path = fig_path, width = 21, height = 15, units = "cm")

corr_sum_plot <- summary_plot(filter(metrics_summary, scenario != scenario_labels[["1"]]),
                              "mean_corr", "mcse_corr",
                              "Mean correlation with the true CATE", "Mean correlation")
ggsave("cts_corr_summary.png", plot = corr_sum_plot, path = fig_path,
       width = 21, height = 15, units = "cm")

# --- BLP test p-values --------------------------------------------------
BLP_plot <- metrics %>%
  ggplot(aes(x = n, y = BLP_p, colour = model)) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_wrap(~scenario, scales = "free_x") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "BLP test p-values", y = "p-value", x = "Sample size", colour = "Model")
ggsave("cts_blp_all.png", path = fig_path, width = 21, height = 15, units = "cm")

BLP_sum_plot <- summary_plot(metrics_summary, "mean_BLP", "mcse_BLP",
                             "Mean BLP test p-value", "Mean p-value", hline = 0.05)
ggsave("cts_blp_summary.png", plot = BLP_sum_plot, path = fig_path,
       width = 21, height = 15, units = "cm")

# --- CATE permutation test p-values --------------------------------------
indep_cate_plot <- metrics %>%
  ggplot(aes(x = n, y = indep_cate, colour = model)) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_wrap(~scenario, scales = "free_x") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "CATE permutation test p-values", y = "p-value", x = "Sample size",
       colour = "Model")
ggsave("cts_indep_cate_all.png", path = fig_path, width = 21, height = 15, units = "cm")

indep_cate_sum_plot <- summary_plot(metrics_summary, "mean_indep_cate", "mcse_indep_cate",
                                    "Mean CATE permutation test p-value", "Mean p-value",
                                    hline = 0.05)
ggsave("cts_indep_cate_summary.png", plot = indep_cate_sum_plot, path = fig_path,
       width = 21, height = 15, units = "cm")

# --- PO permutation test p-values ----------------------------------------
indep_po_plot <- metrics %>%
  ggplot(aes(x = n, y = indep_po, colour = model)) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_wrap(~scenario, scales = "free_x") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "PO permutation test p-values", y = "p-value", x = "Sample size",
       colour = "Model")
ggsave("cts_indep_po_all.png", path = fig_path, width = 21, height = 15, units = "cm")

indep_po_sum_plot <- summary_plot(metrics_summary, "mean_indep_po", "mcse_indep_po",
                                  "Mean PO permutation test p-value", "Mean p-value",
                                  hline = 0.05)
ggsave("cts_indep_po_summary.png", plot = indep_po_sum_plot, path = fig_path,
       width = 21, height = 15, units = "cm")

# --- headline table -----------------------------------------------------
headline <- metrics_summary %>%
  select(scenario, model, n, mean_bias, mean_mse, mean_corr, mean_BLP,
         mean_indep_cate, mean_indep_po) %>%
  arrange(scenario, n, mean_mse)

print(headline, n = Inf)
saveRDS(headline, file.path(res_path, "cts_headline.RDS"))
write_csv(headline, file.path(res_path, "cts_headline.csv"))

print("figures written!")
