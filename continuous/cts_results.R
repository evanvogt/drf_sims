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
fig_path <- file.path(dirname(path), "results", "all_figures", "continuous")
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
# mcse denominator is sqrt(sum(!is.na(x))), matching the non-NA count of that
# column - not sqrt(n()), the whole group's row count regardless of NAs.
# BLP_p/indep_cate/indep_po are NA_real_ whenever the HTE test wasn't run
# (R/metrics.R::hte_test_metrics()), so those columns are the most exposed.
metrics_summary <- metrics %>%
  group_by(scenario, model, n) %>%
  summarise(
    mean_bias = mean(bias, na.rm = T),
    mcse_bias = sd(bias, na.rm = T) / sqrt(sum(!is.na(bias))),
    mean_ate_bias = mean(ate_bias, na.rm = T),
    mcse_ate_bias = sd(ate_bias, na.rm = T) / sqrt(sum(!is.na(ate_bias))),
    mean_rel_ate_bias = mean(rel_ate_bias, na.rm = T),
    mcse_rel_ate_bias = sd(rel_ate_bias, na.rm = T) / sqrt(sum(!is.na(rel_ate_bias))),
    mean_rel_bias_cate = mean(rel_bias_cate, na.rm = T),
    mcse_rel_bias_cate = sd(rel_bias_cate, na.rm = T) / sqrt(sum(!is.na(rel_bias_cate))),
    mean_mse = mean(mse, na.rm = T),
    mcse_mse = sd(mse, na.rm = T) / sqrt(sum(!is.na(mse))),
    mean_rmse = mean(rmse, na.rm = T),
    mcse_rmse = sd(rmse, na.rm = T) / sqrt(sum(!is.na(rmse))),
    mean_mae = mean(mae, na.rm = T),
    mcse_mae = sd(mae, na.rm = T) / sqrt(sum(!is.na(mae))),
    mean_corr = mean(corr, na.rm = T),
    mcse_corr = sd(corr, na.rm = T) / sqrt(sum(!is.na(corr))),
    mean_spearman = mean(spearman, na.rm = T),
    mcse_spearman = sd(spearman, na.rm = T) / sqrt(sum(!is.na(spearman))),
    mean_sign_acc = mean(sign_acc, na.rm = T),
    mcse_sign_acc = sd(sign_acc, na.rm = T) / sqrt(sum(!is.na(sign_acc))),
    mean_BLP = mean(BLP_p, na.rm = T),
    mcse_BLP = sd(BLP_p, na.rm = T) / sqrt(sum(!is.na(BLP_p))),
    mean_indep_cate = mean(indep_cate, na.rm = T),
    mcse_indep_cate = sd(indep_cate, na.rm = T) / sqrt(sum(!is.na(indep_cate))),
    mean_indep_po = mean(indep_po, na.rm = T),
    mcse_indep_po = sd(indep_po, na.rm = T) / sqrt(sum(!is.na(indep_po))),
    # power: proportion of runs rejecting at alpha=0.05 - a genuine 0/1
    # indicator, so it gets the binomial MCSE, unlike the mean p-values above
    power_BLP = mean(BLP_p < 0.05, na.rm = T),
    mcse_power_BLP = sqrt(power_BLP * (1 - power_BLP) / sum(!is.na(BLP_p))),
    power_indep_cate = mean(indep_cate < 0.05, na.rm = T),
    mcse_power_indep_cate = sqrt(power_indep_cate * (1 - power_indep_cate) / sum(!is.na(indep_cate))),
    power_indep_po = mean(indep_po < 0.05, na.rm = T),
    mcse_power_indep_po = sqrt(power_indep_po * (1 - power_indep_po) / sum(!is.na(indep_po))),
    mean_n_na = mean(n_na, na.rm = T),
    total_n_na = sum(n_na, na.rm = T),
    .groups = "drop"
  )

# helper: the house summary plot, points with a 95% CI (mean +/- 1.96 x MCSE)
# error bar, faceted by scenario (there's no family/variant dimension here -
# just scenario x model x n) - see crossfitting/cf_results.R's summary_plot
# for the pattern this follows.
#
# NOTE: this used to draw +/- 1x MCSE (~68% coverage), which made a
# genuinely unbiased estimator look biased on ~32% of points by chance alone
# (verified against dr_oracle's fresh 100-run continuous results: 42% of
# cells missed zero at +/- 1x MCSE vs 5% - the expected false-positive rate -
# at the proper 95% CI below). qnorm(1 - alpha/2), not a hardcoded 1.96, to
# match R/metrics.R::normal_interval's convention.
summary_plot <- function(df, mean_col, mcse_col, title, ylab, hline = 0, alpha = 0.05) {
  z <- qnorm(1 - alpha / 2)
  df %>%
    mutate(est = .data[[mean_col]],
           lo = .data[[mean_col]] - z * .data[[mcse_col]],
           hi = .data[[mean_col]] + z * .data[[mcse_col]]) %>%
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

# --- ATE bias --------------------------------------------------------------
# bias of the average effect (mean(est) - mean(true)), as opposed to bias's
# per-unit average - can differ in sign/magnitude when errors don't average
# out evenly across units
ate_bias_plot <- metrics %>%
  ggplot(aes(x = n, y = ate_bias, colour = model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_wrap(~scenario, scales = "free") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Bias in the average treatment effect", y = "ATE bias",
       x = "Sample size", colour = "Model")
ggsave("cts_ate_bias_all.png", path = fig_path, width = 21, height = 15, units = "cm")

ate_bias_sum_plot <- summary_plot(metrics_summary, "mean_ate_bias", "mcse_ate_bias",
                                  "Mean ATE bias", "Mean ATE bias")
ggsave("cts_ate_bias_summary.png", plot = ate_bias_sum_plot, path = fig_path,
       width = 21, height = 15, units = "cm")

# --- relative bias -----------------------------------------------------
# rel_ate_bias: ATE bias divided by the true ATE. rel_bias_cate: per-unit
# (est - true) / true, averaged over units (NA wherever true CATE is exactly
# 0 for a unit - see R/metrics.R::cate_metrics())
rel_ate_bias_sum_plot <- summary_plot(metrics_summary, "mean_rel_ate_bias",
                                      "mcse_rel_ate_bias",
                                      "Mean relative bias in the ATE",
                                      "Mean relative ATE bias")
ggsave("cts_rel_ate_bias_summary.png", plot = rel_ate_bias_sum_plot, path = fig_path,
       width = 21, height = 15, units = "cm")

rel_bias_cate_sum_plot <- summary_plot(metrics_summary, "mean_rel_bias_cate",
                                       "mcse_rel_bias_cate",
                                       "Mean relative bias in the CATE",
                                       "Mean relative CATE bias")
ggsave("cts_rel_bias_cate_summary.png", plot = rel_bias_cate_sum_plot, path = fig_path,
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

# --- RMSE ------------------------------------------------------------------
# sqrt(mse) - same ranking as MSE, shown on the outcome's own scale
rmse_plot <- metrics %>%
  ggplot(aes(x = n, y = rmse, colour = model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_wrap(~scenario, scales = "free") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "RMSE of CATE estimates", y = "RMSE", x = "Sample size", colour = "Model")
ggsave("cts_rmse_all.png", path = fig_path, width = 21, height = 15, units = "cm")

rmse_sum_plot <- summary_plot(metrics_summary, "mean_rmse", "mcse_rmse",
                              "Mean RMSE of CATE estimates", "Mean RMSE")
ggsave("cts_rmse_summary.png", plot = rmse_sum_plot, path = fig_path,
       width = 21, height = 15, units = "cm")

# --- MAE -------------------------------------------------------------------
# mean absolute error - a robustness check on MSE, less sensitive to outliers
mae_plot <- metrics %>%
  ggplot(aes(x = n, y = mae, colour = model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_wrap(~scenario, scales = "free") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "MAE of CATE estimates", y = "MAE", x = "Sample size", colour = "Model")
ggsave("cts_mae_all.png", path = fig_path, width = 21, height = 15, units = "cm")

mae_sum_plot <- summary_plot(metrics_summary, "mean_mae", "mcse_mae",
                             "Mean MAE of CATE estimates", "Mean MAE")
ggsave("cts_mae_summary.png", plot = mae_sum_plot, path = fig_path,
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

# --- Spearman rank correlation with the truth -------------------------------
# same scenario-1-undefined convention as Pearson corr (R/metrics.R)
spearman_plot <- corr_df %>%
  ggplot(aes(x = n, y = spearman, colour = model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_wrap(~scenario, scales = "free_y") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Spearman rank correlation between estimated and true CATEs",
       y = "Spearman correlation", x = "Sample size", colour = "Model")
ggsave("cts_spearman_all.png", path = fig_path, width = 21, height = 15, units = "cm")

spearman_sum_plot <- summary_plot(filter(metrics_summary, scenario != scenario_labels[["1"]]),
                                  "mean_spearman", "mcse_spearman",
                                  "Mean Spearman correlation with the true CATE",
                                  "Mean Spearman correlation")
ggsave("cts_spearman_summary.png", plot = spearman_sum_plot, path = fig_path,
       width = 21, height = 15, units = "cm")

# --- sign accuracy -----------------------------------------------------------
# proportion of units where sign(est) == sign(true) - direction detection,
# not magnitude. 0.5 is the chance level for a coin-flip on sign
sign_acc_plot <- metrics %>%
  ggplot(aes(x = n, y = sign_acc, colour = model)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_wrap(~scenario, scales = "free_x") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Sign accuracy of CATE estimates", y = "Proportion correct sign",
       x = "Sample size", colour = "Model")
ggsave("cts_sign_acc_all.png", path = fig_path, width = 21, height = 15, units = "cm")

sign_acc_sum_plot <- summary_plot(metrics_summary, "mean_sign_acc", "mcse_sign_acc",
                                  "Mean sign accuracy", "Mean proportion correct sign",
                                  hline = 0.5)
ggsave("cts_sign_acc_summary.png", plot = sign_acc_sum_plot, path = fig_path,
       width = 21, height = 15, units = "cm")

# --- true CATE test evaluation -----------------------------------------
# BLP and independence tests run on the true CATE and true nuisances instead
# of an estimator's (see continuous/README.md's "True-CATE HTE test
# evaluation") - the ceiling each model's own BLP/indep_cate test is chasing.
# No indep_po counterpart: dr_oracle's own indep_po already covers that (it's
# the true pseudo-outcome test, and dr_oracle is already one of the 5 models
# below). Bound in as an extra "model" only where it's plotted (the BLP/
# indep_cate panels that follow) - every other plot above/below is unaffected.
true_cate_tests <- readRDS(file.path(res_path, "cts_true_cate_tests.RDS")) %>%
  mutate(
    scenario = factor(scenario, levels = names(scenario_labels) %>% as.integer(),
                      labels = unname(scenario_labels)),
    n = factor(n, levels = c(100, 250, 500, 1000))
  ) %>%
  droplevels()

true_cate_summary <- true_cate_tests %>%
  group_by(scenario, n) %>%
  summarise(
    mean_BLP = mean(BLP_p, na.rm = T),
    mcse_BLP = sd(BLP_p, na.rm = T) / sqrt(sum(!is.na(BLP_p))),
    mean_indep_cate = mean(indep_cate, na.rm = T),
    mcse_indep_cate = sd(indep_cate, na.rm = T) / sqrt(sum(!is.na(indep_cate))),
    power_BLP = mean(BLP_p < 0.05, na.rm = T),
    mcse_power_BLP = sqrt(power_BLP * (1 - power_BLP) / sum(!is.na(BLP_p))),
    power_indep_cate = mean(indep_cate < 0.05, na.rm = T),
    mcse_power_indep_cate = sqrt(power_indep_cate * (1 - power_indep_cate) / sum(!is.na(indep_cate))),
    .groups = "drop"
  ) %>%
  mutate(model = "True CATE")

true_cate_raw <- true_cate_tests %>% mutate(model = "True CATE")

model_levels_true <- c(levels(metrics$model), "True CATE")

# row-bind the true-CATE reference series onto a per-run or summarised data
# frame, extending `model`'s levels so it becomes one more colour/facet entry
bind_true_cate <- function(df, true_df, extra_cols) {
  bind_rows(
    mutate(df, model = factor(as.character(model), levels = model_levels_true)),
    true_df %>%
      mutate(model = factor(model, levels = model_levels_true)) %>%
      select(scenario, n, model, all_of(extra_cols))
  )
}

# --- BLP test p-values --------------------------------------------------
BLP_plot <- bind_true_cate(metrics, true_cate_raw, "BLP_p") %>%
  ggplot(aes(x = n, y = BLP_p, colour = model)) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_wrap(~scenario, scales = "free_x") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "BLP test p-values", y = "p-value", x = "Sample size", colour = "Model")
ggsave("cts_blp_all.png", path = fig_path, width = 21, height = 15, units = "cm")

BLP_sum_plot <- summary_plot(
  bind_true_cate(metrics_summary, true_cate_summary, c("mean_BLP", "mcse_BLP")),
  "mean_BLP", "mcse_BLP", "Mean BLP test p-value", "Mean p-value", hline = 0.05)
ggsave("cts_blp_summary.png", plot = BLP_sum_plot, path = fig_path,
       width = 21, height = 15, units = "cm")

# --- CATE permutation test p-values --------------------------------------
indep_cate_plot <- bind_true_cate(metrics, true_cate_raw, "indep_cate") %>%
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

indep_cate_sum_plot <- summary_plot(
  bind_true_cate(metrics_summary, true_cate_summary, c("mean_indep_cate", "mcse_indep_cate")),
  "mean_indep_cate", "mcse_indep_cate",
  "Mean CATE permutation test p-value", "Mean p-value", hline = 0.05)
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

# --- HTE test power (proportion rejecting at alpha=0.05) --------------------
# the rsimsum-standard measure for a hypothesis test's simulated behaviour,
# alongside the mean p-values above
BLP_power_plot <- summary_plot(
  bind_true_cate(metrics_summary, true_cate_summary, c("power_BLP", "mcse_power_BLP")),
  "power_BLP", "mcse_power_BLP", "BLP test power", "Power", hline = 0.05)
ggsave("cts_blp_power.png", plot = BLP_power_plot, path = fig_path,
       width = 21, height = 15, units = "cm")

indep_cate_power_plot <- summary_plot(
  bind_true_cate(metrics_summary, true_cate_summary,
                c("power_indep_cate", "mcse_power_indep_cate")),
  "power_indep_cate", "mcse_power_indep_cate",
  "CATE permutation test power", "Power", hline = 0.05)
ggsave("cts_indep_cate_power.png", plot = indep_cate_power_plot, path = fig_path,
       width = 21, height = 15, units = "cm")

indep_po_power_plot <- summary_plot(metrics_summary, "power_indep_po",
                                    "mcse_power_indep_po",
                                    "PO permutation test power", "Power",
                                    hline = 0.05)
ggsave("cts_indep_po_power.png", plot = indep_po_power_plot, path = fig_path,
       width = 21, height = 15, units = "cm")

# --- missing estimates (n_na) diagnostic ------------------------------------
# a count, mostly zero - not a performance metric, so a table rather than a
# boxplot; printed only, not saved, same treatment as every other
# intermediate table here
na_table <- metrics_summary %>%
  select(scenario, model, n, mean_n_na, total_n_na) %>%
  filter(total_n_na > 0) %>%
  arrange(scenario, n, model)

if (nrow(na_table) > 0) {
  print(na_table, n = Inf)
} else {
  print("no missing CATE estimates in any run")
}

# --- headline table -----------------------------------------------------
headline <- metrics_summary %>%
  select(scenario, model, n, mean_bias, mean_ate_bias, mean_rel_ate_bias,
         mean_rel_bias_cate, mean_mse, mean_rmse, mean_mae,
         mean_sign_acc, mean_corr, mean_BLP, mean_indep_cate, mean_indep_po,
         power_BLP, power_indep_cate, power_indep_po) %>%
  arrange(scenario, n, mean_mse)

print(headline, n = Inf)
saveRDS(headline, file.path(res_path, "cts_headline.RDS"))
write_csv(headline, file.path(res_path, "cts_headline.csv"))

print("figures written!")
