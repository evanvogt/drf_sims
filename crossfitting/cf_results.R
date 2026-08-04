##########
# title: figures for the crossfitting comparison
##########

# libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(paletteer)
library(here)
library(patchwork)
library(purrr)

# functions
source(here("crossfitting", "cf_metrics.R"))

# paths
path <- here()
res_path <- file.path(dirname(path), "results", "crossfitting")
fig_path <- file.path(dirname(path), "results", "thesis_figures", "crossfitting")
dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)

metrics <- readRDS(file.path(res_path, "cf_metrics.RDS"))

# tidy up
metrics <- metrics %>%
  mutate(
    scenario = factor(
      case_when(scenario == 1 ~ "Null",
                scenario == 4 ~ "Simple",
                scenario == 6 ~ "Interaction",
                scenario == 9 ~ "Non-linear"),
      levels = c("Null", "Simple", "Interaction", "Non-linear")),
    family = factor(recode(family, !!!family_labels), levels = unname(family_labels)),
    variant = factor(recode(variant, !!!variant_labels),
                     levels = unname(variant_labels[variant_levels])),
    set = factor(set, levels = c("train", "test"), labels = c("Training sample", "Test sample"))
  ) %>%
  droplevels()

# per variant summaries
metrics_summary <- metrics %>%
  group_by(scenario, family, variant, set) %>%
  summarise(
    mean_bias = mean(bias, na.rm = T),
    mcse_bias = sd(bias, na.rm = T) / sqrt(n()),
    mean_mse = mean(mse, na.rm = T),
    mcse_mse = sd(mse, na.rm = T) / sqrt(n()),
    mean_corr = mean(corr, na.rm = T),
    mcse_corr = sd(corr, na.rm = T) / sqrt(n()),
    mean_sign_acc = mean(sign_acc, na.rm = T),
    mcse_sign_acc = sd(sign_acc, na.rm = T) / sqrt(n()),
    mean_mse_single = mean(mse_test_single, na.rm = T),
    mcse_mse_single = sd(mse_test_single, na.rm = T) / sqrt(n()),
    mean_time = mean(time_total, na.rm = T),
    mean_time_nuisance = mean(time_nuisance, na.rm = T),
    mean_time_stage2 = mean(time_stage2, na.rm = T),
    .groups = "drop"
  )

# helper: the house summary plot, points with MCSE error bars
summary_plot <- function(df, mean_col, mcse_col, title, ylab, hline = 0) {
  df %>%
    mutate(est = .data[[mean_col]],
           lo = .data[[mean_col]] - .data[[mcse_col]],
           hi = .data[[mean_col]] + .data[[mcse_col]]) %>%
    ggplot(aes(x = variant, y = est, colour = set, ymin = lo, ymax = hi)) +
    geom_hline(yintercept = hline, linetype = "dashed") +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbar(position = position_dodge(width = 0.5), linewidth = 0.3, width = 0.3) +
    facet_grid(rows = vars(family), cols = vars(scenario), scales = "free_y") +
    scale_colour_paletteer_d("rcartocolor::Safe") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title, y = ylab, x = "Crossfitting procedure", colour = "Scored on")
}

# --- MSE: the headline ------------------------------------------------------
mse_plot <- metrics %>%
  ggplot(aes(x = variant, y = mse, colour = set)) +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_grid(rows = vars(family), cols = vars(scenario), scales = "free_y") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "MSE of CATE estimates by crossfitting procedure",
       y = "MSE", x = "Crossfitting procedure", colour = "Scored on")
ggsave("cf_mse_all.png", path = fig_path, width = 21, height = 15, units = "cm")

mse_sum_plot <- summary_plot(metrics_summary, "mean_mse", "mcse_mse",
                             "Mean MSE of CATE estimates", "Mean MSE")
ggsave("cf_mse_summary.png", plot = mse_sum_plot, path = fig_path,
       width = 21, height = 15, units = "cm")

corr_sum_plot <- summary_plot(filter(metrics_summary, scenario != "Null"),
                              "mean_corr", "mcse_corr",
                              "Mean correlation with the true CATE", "Mean correlation")
ggsave("cf_corr_summary.png", plot = corr_sum_plot, path = fig_path,
       width = 21, height = 15, units = "cm")

# --- generalisation gap: test MSE minus train MSE, per run ------------------
# NOT an optimism measure. Every arm is scored against the known true CATE, not
# against the pseudo-outcomes it was fit to, so there is no optimism to detect -
# an in-sample prediction is contaminated by the unit's own noisy pseudo-outcome
# and so scores worse against truth, not better. This plot shows how far each
# fitted surface degrades on fresh covariate draws.
gap <- metrics %>%
  select(scenario, family, variant, run, set, mse) %>%
  pivot_wider(names_from = set, values_from = mse) %>%
  mutate(gap = `Test sample` - `Training sample`)

gap_plot <- gap %>%
  ggplot(aes(x = variant, y = gap, colour = family)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_grid(rows = vars(family), cols = vars(scenario), scales = "free") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(title = "Generalisation gap: test MSE minus training MSE",
       subtitle = "How far each fitted surface degrades on fresh covariate draws",
       y = "Test MSE - training MSE", x = "Crossfitting procedure")
ggsave("cf_generalisation_gap.png", path = fig_path, width = 21, height = 15, units = "cm")

# --- ensemble vs single model on the test set -------------------------------
# a crossfit arm deploys V fold models and averages their test predictions, while
# a whole-sample arm has one model. mean_mse is the ensemble; mean_mse_single
# scores each fold model separately, which is the like-for-like reading. the two
# coincide for the whole-sample arms, so any gap is the ensembling effect.
ensemble_plot <- metrics_summary %>%
  filter(set == "Test sample") %>%
  select(scenario, family, variant, mean_mse, mean_mse_single) %>%
  pivot_longer(c(mean_mse, mean_mse_single),
               names_to = "scoring", values_to = "mse") %>%
  mutate(scoring = recode(scoring,
                          mean_mse = "Averaged fold models (deployed)",
                          mean_mse_single = "Single fold model (like-for-like)")) %>%
  ggplot(aes(x = variant, y = mse, colour = scoring)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  facet_grid(rows = vars(family), cols = vars(scenario), scales = "free_y") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Test MSE: ensemble of fold models against a single fold model",
       subtitle = "Identical for whole-sample arms; the gap on crossfit arms is the ensembling effect",
       y = "Mean test MSE", x = "Crossfitting procedure", colour = "Scoring")
ggsave("cf_ensemble_effect.png", path = fig_path, width = 21, height = 15, units = "cm")

# --- bias -------------------------------------------------------------------
bias_plot <- metrics %>%
  ggplot(aes(x = variant, y = bias, colour = set)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_grid(rows = vars(family), cols = vars(scenario), scales = "free_y") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Bias of CATE estimates by crossfitting procedure",
       y = "Average bias (estimate - truth)", x = "Crossfitting procedure",
       colour = "Scored on")
ggsave("cf_bias_all.png", path = fig_path, width = 21, height = 15, units = "cm")

# --- correlation with the truth --------------------------------------------
corr_plot <- metrics %>%
  filter(scenario != "Null") %>%
  ggplot(aes(x = variant, y = corr, colour = set)) +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_grid(rows = vars(family), cols = vars(scenario)) +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation between estimated and true CATEs",
       y = "Pearson correlation", x = "Crossfitting procedure", colour = "Scored on")
ggsave("cf_corr_all.png", path = fig_path, width = 21, height = 15, units = "cm")

# --- cost -------------------------------------------------------------------
# what the extra fold pairs actually buy, in seconds
time_plot <- metrics_summary %>%
  filter(set == "Test sample") %>%
  select(scenario, family, variant, mean_time_nuisance, mean_time_stage2) %>%
  pivot_longer(c(mean_time_nuisance, mean_time_stage2),
               names_to = "stage", values_to = "seconds") %>%
  mutate(stage = recode(stage,
                        mean_time_nuisance = "Stage 1 (nuisances)",
                        mean_time_stage2 = "Stage 2 (final model)")) %>%
  ggplot(aes(x = variant, y = seconds, fill = stage)) +
  geom_col() +
  facet_grid(rows = vars(family), cols = vars(scenario), scales = "free_y") +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Mean runtime per replicate",
       y = "Seconds", x = "Crossfitting procedure", fill = "Stage")
ggsave("cf_runtime.png", path = fig_path, width = 21, height = 15, units = "cm")

# --- accuracy against cost --------------------------------------------------
# the decision plot: is double crossfitting on the efficient frontier?
frontier_plot <- metrics_summary %>%
  filter(set == "Test sample") %>%
  ggplot(aes(x = mean_time, y = mean_mse, colour = variant, shape = family)) +
  geom_point(size = 3) +
  facet_wrap(~scenario, scales = "free") +
  scale_x_log10() +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(title = "Test-sample MSE against runtime",
       subtitle = "Down and to the left is better",
       y = "Mean test MSE", x = "Mean runtime per replicate (s, log scale)",
       colour = "Procedure", shape = "Estimator")
ggsave("cf_frontier.png", path = fig_path, width = 21, height = 15, units = "cm")

# --- headline table ---------------------------------------------------------
headline <- metrics_summary %>%
  filter(set == "Test sample") %>%
  select(scenario, family, variant, mean_mse, mcse_mse, mean_mse_single,
         mean_corr, mean_time) %>%
  arrange(family, scenario, mean_mse)

print(headline, n = Inf)
saveRDS(headline, file.path(res_path, "cf_headline.RDS"))

print("figures written!")
