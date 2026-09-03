################
# Title: ICTMC Figures
#################

library(here)
library(tidyverse)
library(paletteer)
library(patchwork)
library(ggridges)

path <- here()
figdir <- file.path(dirname(path), "results", "ICTMC_figs")
dir.create(figdir, showWarnings = F)

metrics_all <- readRDS(file.path(
  dirname(path),
  "collected_metrics",
  "all_metrics.RDS"
))


# Labels and things ----------
scenario_labels <- c(
  `1` = "No HTE",
  `2` = "Binary HTE",
  `3` = "Continuous HTE",
  `4` = "Two HTE vars",
  `5` = "Cts x binary interaction",
  `6` = "Two vars + cts x binary",
  `7` = "Cts x cts interaction",
  `8` = "Two vars + cts x cts",
  `9` = "Cosine HTE",
  `10` = "Exponential HTE"
)

model_labels <- c(
  causal_forest = "Causal forest",
  dr_random_forest = "DR-RandomForest",
  dr_semi_oracle = "DR-semi-oracle",
  dr_superlearner = "DR-SuperLearner"
)

test_labels <- c(
  BLP_p = "BLP",
  indep_po = "Permutation PO",
  indep_cate = "Permutation CATE"
)

# colour palette
model_palette <- c(
  "Causal forest" = "#008080",
  "DR-RandomForest" = "#c71585",
  "DR-semi-oracle" = "#4b0082",
  "DR-SuperLearner" = "#ff4500"
)

test_palette <- c(
  "BLP" = "#008080",
  "Permutation PO" = "#4b0082",
  "Permutation CATE" = "#ff4500"
)

# Sample size and HTE testing ----------

metrics_ss <- metrics_all %>%
  filter(study_name %in% c("binary", "continuous")) %>%
  filter(model != "dr_oracle") %>%
  # keep scenarios 1, 3, 8, 9
  filter(scenario %in% c(1, 3, 8, 9)) %>%
  select(
    scenario,
    n,
    run,
    model,
    bias,
    ate_bias,
    mse,
    rmse,
    mae,
    corr,
    spearman,
    sign_acc,
    n_na,
    BLP_p,
    indep_cate,
    indep_po,
    study_name,
    category
  ) %>%
  mutate(
    scenario = factor(
      scenario,
      levels = names(scenario_labels) %>% as.integer(),
      labels = unname(scenario_labels)
    ),
    n = factor(n, levels = c(100, 250, 500, 1000)),
    model = factor(
      recode(model, !!!model_labels),
      levels = unname(model_labels)
    )
  ) %>%
  droplevels()

met_sum_ss <- metrics_ss %>%
  group_by(scenario, n, model, study_name) %>%
  summarise(
    mean_bias = mean(bias, na.rm = T),
    mcse_bias = sd(bias, na.rm = T) / sqrt(n()),
    mean_ate_bias = mean(ate_bias, na.rm = T),
    mcse_ate_bias = sd(ate_bias, na.rm = T) / sqrt(n()),
    mean_mse = mean(mse, na.rm = T),
    mcse_mse = sd(mse, na.rm = T) / sqrt(n()),
    mean_rmse = mean(rmse, na.rm = T),
    mcse_rmse = sd(rmse, na.rm = T) / sqrt(n()),
    mean_mae = mean(mae, na.rm = T),
    mcse_mae = sd(mae, na.rm = T) / sqrt(n()),
    mean_corr = mean(corr, na.rm = T),
    mcse_corr = sd(corr, na.rm = T) / sqrt(n()),
    mean_spearman = mean(spearman, na.rm = T),
    mcse_spearman = sd(spearman, na.rm = T) / sqrt(n()),
    mean_sign_acc = mean(sign_acc, na.rm = T),
    mcse_sign_acc = sd(sign_acc, na.rm = T) / sqrt(n()),
    mean_BLP = mean(BLP_p, na.rm = T),
    mcse_BLP = sd(BLP_p, na.rm = T) / sqrt(n()),
    mean_indep_cate = mean(indep_cate, na.rm = T),
    mcse_indep_cate = sd(indep_cate, na.rm = T) / sqrt(n()),
    mean_indep_po = mean(indep_po, na.rm = T),
    mcse_indep_po = sd(indep_po, na.rm = T) / sqrt(n()),
    mean_n_na = mean(n_na, na.rm = T),
    total_n_na = sum(n_na, na.rm = T),
    .groups = "drop"
  )


summary_plot <- function(df, mean_col, mcse_col, title, ylab, hline = 0) {
  df %>%
    mutate(
      est = .data[[mean_col]],
      lo = .data[[mean_col]] - .data[[mcse_col]],
      hi = .data[[mean_col]] + .data[[mcse_col]]
    ) %>%
    ggplot(aes(x = n, y = est, colour = model, ymin = lo, ymax = hi)) +
    geom_hline(yintercept = hline, linetype = "dashed") +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbar(
      position = position_dodge(width = 0.5),
      linewidth = 0.3,
      width = 0.3
    ) +
    facet_grid(
      cols = vars(scenario),
      rows = vars(study_name),
      scales = "free"
    ) +
    scale_colour_manual(values = model_palette) +
    theme_light() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(colour = "black"),
      legend.position = "bottom"
    ) +
    labs(title = title, y = ylab, x = "Sample size", colour = "Model")
}


bias_ss_plot <- summary_plot(
  met_sum_ss,
  "mean_bias",
  "mcse_bias",
  "Mean Bias in CATE Estimates",
  "Mean Bias"
)
ggsave(
  "ss_bias.png",
  bias_ss_plot,
  path = figdir,
  height = 20,
  width = 30,
  units = "cm"
)


mae_ss_plot <- summary_plot(
  met_sum_ss,
  "mean_mae",
  "mcse_mae",
  "Mean Absolute Error in CATE Estimates",
  "MAE"
)

ggsave(
  "ss_mae.png",
  mae_ss_plot,
  path = figdir,
  height = 20,
  width = 30,
  units = "cm"
)


# HTE tests ----------
met_tests <- metrics_ss %>%
  select(scenario, n, run, model, BLP_p, indep_po, indep_cate, study_name) %>%
  pivot_longer(
    c(BLP_p, indep_po, indep_cate),
    names_to = "test",
    values_to = "pval"
  ) %>%
  mutate(
    test = factor(recode(test, !!!test_labels), levels = unname(test_labels))
  )


ridge_plot <- function(df, title) {
  df %>%
    ggplot(aes(x = pval, y = test, colour = test, fill = test)) +
    stat_density_ridges(alpha = 0.5, na.rm = TRUE, from = 0, to = 1) +
    geom_vline(xintercept = 0.1, linetype = "dashed") +
    facet_grid(cols = vars(scenario), rows = vars(n)) +
    scale_colour_manual(values = test_palette) +
    scale_fill_manual(values = test_palette) +
    theme_light() +
    theme(
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(colour = "black"),
      legend.position = "bottom"
    ) +
    labs(title = title, y = NULL) +
    xlim(0, 1)
}

binary_ridge <- ridge_plot(
  filter(met_tests, study_name == "binary"),
  "Binary outcome"
)
continuous_ridge <- ridge_plot(
  filter(met_tests, study_name == "continuous"),
  "Continuous outcome"
)

binary_ridge /
  continuous_ridge +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


met_tests_sum <- met_tests %>%
  group_by(scenario, n, study_name, test) %>%
  summarise(
    runs = sum(!is.na(pval)),
    mean_p = mean(pval, na.rm = T),
    mcse_p = sqrt(mean_p * (1 - mean_p) / runs),
    .groups = "drop"
  )

met_tests_sum %>%
  mutate(
    hi = mean_p + mcse_p,
    lo = mean_p - mcse_p
  ) %>%
  ggplot(aes(x = n, y = mean_p, ymin = lo, ymax = hi, colour = test)) +
  geom_hline(yintercept = 0.1, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(
    position = position_dodge(width = 0.5),
    linewidth = 0.3,
    width = 0.3
  ) +
  facet_grid(rows = vars(study_name), cols = vars(scenario)) +
  scale_colour_manual(values = test_palette) +
  theme_light() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(colour = "black"),
    legend.position = "bottom"
  ) +
  labs(
    title = "HTE testing p-values",
    y = "p",
    x = "Sample size",
    colour = "Test"
  )
ggsave("test_p.png", path = figdir, height = 20, width = 30, units = "cm")


# missing data -----------------------

metrics_miss <- metrics_all %>%
  filter(study_name %in% c("missing/binary", "missing/continuous")) %>%
  #ilter(scenario %in% c(1, 3, 8, 9)) %>%
  select(
    scenario,
    run,
    model,
    bias,
    ate_bias,
    mse,
    rmse,
    mae,
    corr,
    spearman,
    sign_acc,
    mechanism,
    method,
    mse_complete,
    rel_efficiency,
    study_name
  ) %>%
  filter(model != "dr_oracle")

met_sum_miss <- metrics_miss %>%
  group_by(scenario, mechanism, study_name, method) %>%
  summarise(
    mean_bias = mean(bias, na.rm = T),
    mcse_bias = sd(bias, na.rm = T) / sqrt(n()),
    mean_ate_bias = mean(ate_bias, na.rm = T),
    mcse_ate_bias = sd(ate_bias, na.rm = T) / sqrt(n()),
    mean_mse = mean(mse, na.rm = T),
    mcse_mse = sd(mse, na.rm = T) / sqrt(n()),
    mean_rmse = mean(rmse, na.rm = T),
    mcse_rmse = sd(rmse, na.rm = T) / sqrt(n()),
    mean_mae = mean(mae, na.rm = T),
    mcse_mae = sd(mae, na.rm = T) / sqrt(n()),
    mean_corr = mean(corr, na.rm = T),
    mcse_corr = sd(corr, na.rm = T) / sqrt(n()),
    mean_spearman = mean(spearman, na.rm = T),
    mcse_spearman = sd(spearman, na.rm = T) / sqrt(n()),
    mean_sign_acc = mean(sign_acc, na.rm = T),
    mcse_sign_acc = sd(sign_acc, na.rm = T) / sqrt(n()),
    mean_rel_eff = mean(rel_efficiency, na.rm = T),
    mcse_rel_eff = sd(rel_efficiency, na.rm = T) / sqrt(n()),
    .groups = "drop"
  )


met_sum_miss %>%
  mutate(
    hi = mean_rel_eff + mcse_rel_eff,
    lo = mean_rel_eff - mcse_rel_eff
  ) %>%
  ggplot(aes(
    x = mechanism,
    y = mean_rel_eff,
    ymin = lo,
    ymax = hi,
    colour = method
  )) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(
    position = position_dodge(width = 0.5),
    linewidth = 0.3,
    width = 0.3
  ) +
  facet_grid(rows = vars(study_name), cols = vars(scenario), scales = "free") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_light() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(colour = "black"),
    legend.position = "bottom"
  ) +
  labs(
    title = "Relative Efficieny of Missing Data Handling Methods",
    y = "Rekative Efficiency",
    x = "Sample size",
    colour = "Method"
  )
