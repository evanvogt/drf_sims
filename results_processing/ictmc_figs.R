################
# Title: ICTMC Figures
#################

library(here)
library(tidyverse)
library(paletteer)
library(patchwork)
library(ggridges)
source(here("R", "figures.R"))

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

# true-value HTE test ceiling panel: BLP/PO computed on the true CATE / true
# pseudo-outcome instead of an estimator's fit (see "True-value HTE tests"
# section below)
true_test_palette <- c(
  "True BLP" = "#4b0082",
  "True PO permutation test" = "#dc143c"
)

# The poster's full colour vocabulary, verbatim from
# ICTMC_poster_planning.md's "expanded palette of an additional 27 colours".
# #0000cd (Imperial blue) plus core black (#161a1d) and white (#ffffff) are
# the poster's chrome (title, logo, background) - reserved, not used for data
# below. model_palette above already happens to be 4 colours drawn from this
# list; the palettes below extend the same approach to every other
# categorical colour dimension the poster's remaining panels need. Colours
# are reused across *different* panels once the 27-colour budget runs thin -
# never within the same legend - matching the doc's "different plots will
# have different numbers of colours".
IMPERIAL_PALETTE <- c(
  "#232333",
  "#000080",
  "#8b4513",
  "#008080",
  "#c71585",
  "#4b0082",
  "#dc143c",
  "#ff4500",
  "#006400",
  "#708090",
  "#0000cd",
  "#ffff00",
  "#40e0d0",
  "#ee82ee",
  "#7b68ee",
  "#ff0000",
  "#ff8c00",
  "#00ff7f",
  "#f5f5f5",
  "#00bfff",
  "#f0e68c",
  "#afeeee",
  "#ffb6c1",
  "#e6e6fa",
  "#fa8072",
  "#ffa500",
  "#98fb98"
)

# missing-data handling method (7 levels, METHOD_LABELS' display names, minus
# single mean/forest imputation which are dropped for the poster - MPA/MIA are
# poster-only abbreviations for "Missing indicators"/"Inbuilt missingness
# handling", applied locally below rather than in the shared METHOD_LABELS)
method_palette <- c(
  "Complete case" = "#708090",
  "Single model-based imputation" = "#7b68ee",
  "MPA" = "#ff8c00",
  "IPW" = "#00bfff",
  "Multiple imputation (rf)" = "#dc143c",
  "MIA" = "#8b4513",
  "Complete data (reference)" = "#232333"
)

# CI panel: two independent dimensions (see cts_ci_results.qmd's `interval`
# derivation, and confidence_intervals/continuous/cts_ci_metrics.R).
# `interval` (colour) - which CI is being scored: the half-sample bootstrap
# (all 4 CI_MODELS) or causal forest's own inbuilt variance estimate
# (causal_forest only, no other model has this row at all).
# `location` (shape) - where it was evaluated: at the sample's own points
# ("Within-sample", every row has this) or at the fixed query grid ("Grid",
# `<model>_grid` rows - exists for all 4 models too, NOT causal_forest-only;
# the inbuilt variance estimate has no grid counterpart, so "CF variance
# estimate" never pairs with "Grid").
ci_interval_palette <- c(
  "Half-sample bootstrap" = "#000080",
  "CF variance estimate" = "#ff0000"
)
ci_location_palette <- c(
  "Within-sample" = 16,
  "Grid" = 17
)

# model evaluation regret panel: nuisance pipeline x score family
proxy_combo_palette <- c(
  "DR risk / AutoML" = "#ffa500",
  "DR risk / XGBoost" = "#40e0d0",
  "Influence / AutoML" = "#00ff7f",
  "Influence / XGBoost" = "#ee82ee"
)

# subgroup validation panel: significance threshold
threshold_palette <- c(
  "0.05" = "#dc143c",
  "0.10" = "#4b0082"
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
    rel_ate_bias,
    rel_bias_cate,
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

# mcse denominator is sqrt(sum(!is.na(x))), matching the non-NA count of that
# column - not sqrt(n()), the whole group's row count regardless of NAs.
# BLP_p/indep_cate/indep_po are NA_real_ whenever the HTE test wasn't run
# (R/metrics.R::hte_test_metrics()), so those columns are the most exposed.
met_sum_ss <- metrics_ss %>%
  group_by(scenario, n, model, study_name) %>%
  summarise(
    mean_bias = mean(bias, na.rm = T),
    mcse_bias = sd(bias, na.rm = T) / sqrt(sum(!is.na(bias))),
    mean_ate_bias = mean(ate_bias, na.rm = T),
    mcse_ate_bias = sd(ate_bias, na.rm = T) / sqrt(sum(!is.na(ate_bias))),
    mean_rel_ate_bias = mean(rel_ate_bias, na.rm = T),
    mcse_rel_ate_bias = sd(rel_ate_bias, na.rm = T) /
      sqrt(sum(!is.na(rel_ate_bias))),
    mean_rel_bias_cate = mean(rel_bias_cate, na.rm = T),
    mcse_rel_bias_cate = sd(rel_bias_cate, na.rm = T) /
      sqrt(sum(!is.na(rel_bias_cate))),
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
    mcse_power_indep_cate = sqrt(
      power_indep_cate * (1 - power_indep_cate) / sum(!is.na(indep_cate))
    ),
    power_indep_po = mean(indep_po < 0.05, na.rm = T),
    mcse_power_indep_po = sqrt(
      power_indep_po * (1 - power_indep_po) / sum(!is.na(indep_po))
    ),
    mean_n_na = mean(n_na, na.rm = T),
    total_n_na = sum(n_na, na.rm = T),
    .groups = "drop"
  )


# 95% CI (mean +/- qnorm(1 - alpha/2) x MCSE) error bar, not a raw +/- 1x
# MCSE (~68% coverage) - see continuous/cts_results.R's summary_plot
summary_plot <- function(
  df,
  mean_col,
  mcse_col,
  title,
  ylab,
  hline = 0,
  alpha = 0.05
) {
  z <- qnorm(1 - alpha / 2)
  df %>%
    mutate(
      est = .data[[mean_col]],
      lo = .data[[mean_col]] - z * .data[[mcse_col]],
      hi = .data[[mean_col]] + z * .data[[mcse_col]]
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

# rel_ate_bias: ATE bias divided by the true ATE. rel_bias_cate: per-unit
# (est - true) / true, averaged over units (NA wherever true CATE is exactly
# 0 for a unit - see R/metrics.R::cate_metrics())
rel_ate_bias_ss_plot <- summary_plot(
  met_sum_ss,
  "mean_rel_ate_bias",
  "mcse_rel_ate_bias",
  "Mean Relative Bias in the ATE",
  "Mean relative ATE bias"
)
ggsave(
  "ss_rel_ate_bias.png",
  rel_ate_bias_ss_plot,
  path = figdir,
  height = 20,
  width = 30,
  units = "cm"
)

rel_bias_cate_ss_plot <- summary_plot(
  met_sum_ss,
  "mean_rel_bias_cate",
  "mcse_rel_bias_cate",
  "Mean Relative Bias in the CATE",
  "Mean relative CATE bias"
)
ggsave(
  "ss_rel_bias_cate.png",
  rel_bias_cate_ss_plot,
  path = figdir,
  height = 20,
  width = 30,
  units = "cm"
)

# Poster panel (final): MAE top / MSE bottom, point+line, faded MCSE CI.
# Continuous outcome only (ICTMC_poster_planning.md's "restrict to continuous
# outcome setting only"), colour keeps the poster's curated model_palette
# rather than point_range_plot()'s thesis-default drf_scale().
met_sum_ss_cts <- filter(met_sum_ss, study_name == "continuous")

ss_mae_panel <- point_range_plot(
  met_sum_ss_cts,
  "mae",
  "Mean absolute error (MAE)",
  x = "n",
  colour = "model",
  facet_rows = NULL,
  facet_cols = "scenario",
  blank_x = TRUE,
  line = TRUE,
  ci_alpha = 0.5,
  palette = scale_colour_manual(values = model_palette)
) +
  labs(x = NULL)

ss_bias_panel <- point_range_plot(
  met_sum_ss_cts,
  "bias",
  "Mean Bias",
  x = "n",
  colour = "model",
  facet_rows = NULL,
  facet_cols = "scenario",
  blank_x = TRUE,
  line = TRUE,
  ci_alpha = 0.5,
  palette = scale_colour_manual(values = model_palette)
) +
  labs(x = NULL)

ss_mse_panel <- point_range_plot(
  met_sum_ss_cts,
  "mse",
  "Mean PEHE (MSE)",
  x = "n",
  colour = "model",
  facet_rows = NULL,
  facet_cols = "scenario",
  facet_scales = "free_y",
  line = TRUE,
  ci_alpha = 0.5,
  palette = scale_colour_manual(values = model_palette)
) +
  labs(x = "Sample size")

(ss_mae_panel / ss_mse_panel) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(
  "ss_mae_mse.png",
  path = figdir,
  height = 20,
  width = 30,
  units = "cm"
)

(ss_bias_panel / ss_mse_panel) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave("ss_bias_mse.png", path = figdir, height = 20, width = 30, units = "cm")


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


# mean_p is the mean of a CONTINUOUS p-value, so its MCSE is the general
# sd(x)/sqrt(n) formula - the binomial p(1-p)/n form (previously used here)
# is only valid for a genuine 0/1 indicator, not a continuous mean, and was
# mathematically wrong. `power` is that genuine 0/1 indicator (reject at
# alpha=0.05), so it correctly gets the binomial formula.
met_tests_sum <- met_tests %>%
  group_by(scenario, n, study_name, test) %>%
  summarise(
    runs = sum(!is.na(pval)),
    mean_p = mean(pval, na.rm = T),
    mcse_p = sd(pval, na.rm = T) / sqrt(runs),
    power_05 = mean(pval < 0.05, na.rm = T),
    mcse_power_05 = sqrt(power_05 * (1 - power_05) / runs),
    power_10 = mean(pval < 0.10, na.rm = T),
    mcse_power_10 = sqrt(power_10 * (1 - power_10) / runs),
    .groups = "drop"
  )

# 95% CI (mean +/- qnorm(0.975) x MCSE), not a raw +/- 1x MCSE (~68% coverage)
met_tests_sum %>%
  mutate(
    hi = mean_p + qnorm(0.975) * mcse_p,
    lo = mean_p - qnorm(0.975) * mcse_p
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

# power: proportion of runs rejecting at alpha=0.05/0.10 - the rsimsum-standard
# measure for a hypothesis test's simulated behaviour, alongside the mean
# p-value above (binomial 95% CI, not the continuous-measure one). Reshaped
# long so threshold becomes a plotted aesthetic (shape) alongside test
# (colour), rather than a fixed alpha baked into the summary.
met_tests_power <- met_tests_sum %>%
  select(
    scenario,
    n,
    study_name,
    test,
    power_05,
    mcse_power_05,
    power_10,
    mcse_power_10
  ) %>%
  pivot_longer(
    cols = c(power_05, mcse_power_05, power_10, mcse_power_10),
    names_to = c(".value", "threshold"),
    names_pattern = "^(mcse_power|power)_(05|10)$"
  ) %>%
  mutate(
    threshold = factor(
      threshold,
      levels = c("05", "10"),
      labels = c("0.05", "0.10")
    ),
    hi = power + qnorm(0.975) * mcse_power,
    lo = power - qnorm(0.975) * mcse_power
  )

# Reference lines mean different things per scenario: in the null scenario
# (no true HTE) they're the nominal size targets a well-calibrated test
# should sit on (0.05/0.10); everywhere else (true HTE present) they're power
# targets (0.90/0.95). Not keyed by study_name, so the same two lines are
# drawn in both facet rows (binary/continuous) for a given scenario column.
null_scenario <- "No HTE"
power_hlines <- bind_rows(
  tibble(scenario = null_scenario, yintercept = c(0.05, 0.10)),
  tibble(
    scenario = setdiff(levels(met_tests_power$scenario), null_scenario),
    yintercept = 0.90
  ),
  tibble(
    scenario = setdiff(levels(met_tests_power$scenario), null_scenario),
    yintercept = 0.95
  )
) %>%
  mutate(scenario = factor(scenario, levels = levels(met_tests_power$scenario)))

met_tests_power %>%
  ggplot(aes(
    x = n,
    y = power,
    ymin = lo,
    ymax = hi,
    colour = test,
    shape = threshold,
    group = interaction(test, threshold)
  )) +
  geom_hline(
    data = power_hlines,
    aes(yintercept = yintercept),
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  geom_line(position = position_dodge(width = 0.6), linewidth = 0.4) +
  geom_point(position = position_dodge(width = 0.6), size = 2) +
  geom_errorbar(
    position = position_dodge(width = 0.6),
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
    title = "HTE testing power (proportion rejecting at alpha=0.05/0.10)",
    y = "Power",
    x = "Sample size",
    colour = "Test",
    shape = "Significance threshold"
  )
ggsave("test_power.png", path = figdir, height = 20, width = 30, units = "cm")


# True-value HTE tests ----------
# Two "ceiling" series, neither tied to any fitted CATE estimator:
#   - True BLP: the BLP test run on the actual ground-truth CATE and true
#     nuisances instead of an estimator's fit (run_true_cate_tests()/
#     true_cate_test_row(), R/cate_models.R). Saved per-study as
#     <prefix>_true_cate_tests.RDS directly in that study's res_path -
#     collect_all_metrics.R's glob only picks up *_metrics.RDS, so this file
#     is never collected into all_metrics.RDS and metrics_all never sees it.
#     Read directly here, same reasoning as val_metrics above (a per-study
#     file that's excluded from the cross-study merge).
#   - True PO permutation test: no true-CATE counterpart needed, since
#     dr_oracle's own indep_po already tests the true pseudo-outcome (known
#     propensity, true outcome model) against X - already an ordinary row of
#     metrics_all, just excluded from metrics_ss's `model != "dr_oracle"`
#     filter above.
read_true_cate_tests <- function(study_dir, prefix, study_label) {
  f <- file.path(
    dirname(path),
    "results",
    study_dir,
    paste0(prefix, "_true_cate_tests.RDS")
  )
  if (!file.exists(f)) {
    message(
      "ictmc_figs: ",
      f,
      " not found locally (not yet synced from the HPC) - ",
      "True BLP will be empty for '",
      study_label,
      "' in this smoke test."
    )
    return(tibble(scenario = integer(), n = double(), BLP_p = double()))
  }
  readRDS(f) %>% mutate(study_name = study_label)
}

true_blp <- bind_rows(
  read_true_cate_tests("binary", "bin", "binary"),
  read_true_cate_tests("continuous", "cts", "continuous")
) %>%
  filter(scenario %in% c(1, 3, 8, 9)) %>%
  transmute(
    scenario = factor(
      scenario,
      levels = names(scenario_labels) %>% as.integer(),
      labels = unname(scenario_labels)
    ),
    n = factor(n, levels = c(100, 250, 500, 1000)),
    study_name,
    test = "True BLP",
    pval = BLP_p
  )

true_po <- metrics_all %>%
  filter(
    study_name %in% c("binary", "continuous"),
    model == "dr_oracle",
    scenario %in% c(1, 3, 8, 9)
  ) %>%
  transmute(
    scenario = factor(
      scenario,
      levels = names(scenario_labels) %>% as.integer(),
      labels = unname(scenario_labels)
    ),
    n = factor(n, levels = c(100, 250, 500, 1000)),
    study_name,
    test = "True PO permutation test",
    pval = indep_po
  )

# no droplevels() - keeping the full 4-level scenario factor even where
# true_blp contributes 0 rows keeps this plot's facet grid aligned with
# test_power.png's
true_tests <- bind_rows(true_blp, true_po) %>%
  mutate(test = factor(test, levels = names(true_test_palette)))

# same power-at-both-thresholds shape as met_tests_power above
true_tests_power <- true_tests %>%
  group_by(scenario, n, study_name, test) %>%
  summarise(
    runs = sum(!is.na(pval)),
    power_05 = mean(pval < 0.05, na.rm = TRUE),
    mcse_power_05 = sqrt(power_05 * (1 - power_05) / runs),
    power_10 = mean(pval < 0.10, na.rm = TRUE),
    mcse_power_10 = sqrt(power_10 * (1 - power_10) / runs),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(power_05, mcse_power_05, power_10, mcse_power_10),
    names_to = c(".value", "threshold"),
    names_pattern = "^(mcse_power|power)_(05|10)$"
  ) %>%
  mutate(
    threshold = factor(
      threshold,
      levels = c("05", "10"),
      labels = c("0.05", "0.10")
    ),
    hi = power + qnorm(0.975) * mcse_power,
    lo = power - qnorm(0.975) * mcse_power
  )

# power_hlines (built above for test_power.png) is keyed only by scenario -
# same nominal-size-vs-power-target meaning applies here unchanged
true_tests_power %>%
  ggplot(aes(
    x = n,
    y = power,
    ymin = lo,
    ymax = hi,
    colour = test,
    shape = threshold,
    group = interaction(test, threshold)
  )) +
  geom_hline(
    data = power_hlines,
    aes(yintercept = yintercept),
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  geom_line(position = position_dodge(width = 0.5), linewidth = 0.4) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(
    position = position_dodge(width = 0.5),
    linewidth = 0.3,
    width = 0.3
  ) +
  facet_grid(rows = vars(study_name), cols = vars(scenario)) +
  scale_colour_manual(values = true_test_palette) +
  theme_light() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(colour = "black"),
    legend.position = "bottom"
  ) +
  labs(
    y = "Power",
    x = "Sample size",
    colour = "Test",
    shape = "Significance threshold"
  )
ggsave(
  "true_test_power.png",
  path = figdir,
  height = 20,
  width = 30,
  units = "cm"
)


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
    rel_ate_bias,
    rel_bias_cate,
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
  filter(model != "dr_oracle") %>%
  # poster only: drop single mean/forest imputation, too cluttered for the
  # panels - thesis chapters keep all 9 methods via the shared METHOD_LABELS
  filter(!method %in% c("mean_imputation", "missforest"))

# unlike metrics_ss, this was never relabelled - apply_labels() (R/figures.R)
# recodes model/method/mechanism/scenario to their display labels in one call
metrics_miss <- metrics_miss %>%
  apply_labels(MISS_SCENARIO_LABELS) %>%
  # poster-only abbreviations to fit panel space - not in the shared
  # METHOD_LABELS, so thesis/qmd output keeps the full names
  mutate(
    method = fct_recode(
      method,
      MPA = "Missing indicators",
      MIA = "Inbuilt missingness handling"
    )
  )

# mcse denominator is sqrt(sum(!is.na(x))), matching the non-NA count of that
# column - not sqrt(n()), the whole group's row count regardless of NAs.
met_sum_miss <- metrics_miss %>%
  group_by(scenario, mechanism, study_name, method) %>%
  summarise(
    mean_bias = mean(bias, na.rm = T),
    mcse_bias = sd(bias, na.rm = T) / sqrt(sum(!is.na(bias))),
    mean_ate_bias = mean(ate_bias, na.rm = T),
    mcse_ate_bias = sd(ate_bias, na.rm = T) / sqrt(sum(!is.na(ate_bias))),
    mean_rel_ate_bias = mean(rel_ate_bias, na.rm = T),
    mcse_rel_ate_bias = sd(rel_ate_bias, na.rm = T) /
      sqrt(sum(!is.na(rel_ate_bias))),
    mean_rel_bias_cate = mean(rel_bias_cate, na.rm = T),
    mcse_rel_bias_cate = sd(rel_bias_cate, na.rm = T) /
      sqrt(sum(!is.na(rel_bias_cate))),
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
    mean_rel_eff = mean(rel_efficiency, na.rm = T),
    mcse_rel_eff = sd(rel_efficiency, na.rm = T) /
      sqrt(sum(!is.na(rel_efficiency))),
    .groups = "drop"
  )


# 95% CI (mean +/- qnorm(0.975) x MCSE), not a raw +/- 1x MCSE (~68% coverage)
met_sum_miss %>%
  mutate(
    hi = mean_rel_eff + qnorm(0.975) * mcse_rel_eff,
    lo = mean_rel_eff - qnorm(0.975) * mcse_rel_eff
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

# Poster panel (final): MAE top / MSE bottom, faceted by scenario, averaged
# over CATE model (met_sum_miss already does this - model isn't in its
# group_by()). Continuous outcome only. Colour is method, for which there's
# no poster-specific palette yet, so this keeps point_range_plot()'s default
# drf_scale() (rcartocolor::Safe - colourblind-friendly, unlike the poster's
# hand-picked 27-colour list per ICTMC_poster_planning.md).
#
# Two variants because whether a connecting line across mechanism
# (MAR/MNAR/MNAR-Y - not a numeric progression, unlike sample size) helps or
# implies a misleading trend is still an open call - see the planning doc.
met_sum_miss_cts <- filter(met_sum_miss, study_name == "missing/continuous")

build_miss_mae_mse_panel <- function(line) {
  top <- point_range_plot(
    met_sum_miss_cts,
    "mae",
    "Mean absolute error (MAE)",
    x = "mechanism",
    colour = "method",
    facet_rows = NULL,
    facet_cols = "scenario",
    # match the bottom (MSE) panel's facet_scales: the default "free" frees
    # both axes per scenario column, and since the "Null" scenario never runs
    # MNAR-Y (R/dgm_scenarios.R), that column then renders a different/
    # narrower x-axis than the others. "free_y" keeps x shared/fixed across
    # columns (only y varies), aligning "Null" with its neighbours and with
    # the MSE panel below.
    facet_scales = "free_y",
    blank_x = TRUE,
    line = line,
    ci_alpha = 0.5,
    palette = scale_colour_manual(values = method_palette)
  ) +
    labs(x = NULL)

  bottom <- point_range_plot(
    met_sum_miss_cts,
    "mse",
    "Mean PEHE (MSE)",
    x = "mechanism",
    colour = "method",
    facet_rows = NULL,
    facet_cols = "scenario",
    facet_scales = "free_y",
    line = line,
    ci_alpha = 0.5,
    palette = scale_colour_manual(values = method_palette)
  ) +
    labs(x = "Missing data mechanism")

  (top / bottom) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

build_miss_mae_mse_panel(line = TRUE)
ggsave(
  "miss_mae_mse_line.png",
  path = figdir,
  height = 20,
  width = 30,
  units = "cm"
)

build_miss_mae_mse_panel(line = FALSE)
ggsave(
  "miss_mae_mse_noline.png",
  path = figdir,
  height = 20,
  width = 30,
  units = "cm"
)

# Same panel shape as build_miss_mae_mse_panel() above, but Bias top / MSE
# bottom instead of MAE top / MSE bottom - hline defaults to 0 in
# point_range_plot(), which is the right reference line for bias.
build_miss_bias_mse_panel <- function(line) {
  top <- point_range_plot(
    met_sum_miss_cts,
    "bias",
    "Bias",
    x = "mechanism",
    colour = "method",
    facet_rows = NULL,
    facet_cols = "scenario",
    # see build_miss_mae_mse_panel()'s comment: keeps x shared/fixed across
    # scenario columns (only y varies) so "Null" (no MNAR-Y) aligns with its
    # neighbours and with the MSE panel below.
    facet_scales = "free_y",
    blank_x = TRUE,
    line = line,
    ci_alpha = 0.5,
    palette = scale_colour_manual(values = method_palette)
  ) +
    labs(x = NULL)

  bottom <- point_range_plot(
    met_sum_miss_cts,
    "mse",
    "Mean PEHE (MSE)",
    x = "mechanism",
    colour = "method",
    facet_rows = NULL,
    facet_cols = "scenario",
    facet_scales = "free_y",
    line = line,
    ci_alpha = 0.5,
    palette = scale_colour_manual(values = method_palette)
  ) +
    labs(x = "Missing data mechanism")

  (top / bottom) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

build_miss_bias_mse_panel(line = TRUE)
ggsave(
  "miss_bias_mse_line.png",
  path = figdir,
  height = 20,
  width = 30,
  units = "cm"
)

build_miss_bias_mse_panel(line = FALSE)
ggsave(
  "miss_bias_mse_noline.png",
  path = figdir,
  height = 20,
  width = 30,
  units = "cm"
)


# Confidence intervals ----------
# confidence_intervals/continuous/cts_ci_metrics.R (via R/metrics.R) writes
# alternative intervals as extra *rows* of `model`: "<model>_inbuilt" for
# causal forest's own pointwise variance (causal_forest only - no other
# model gets this row), "<model>_grid" for the same half-sample bootstrap
# evaluated on the query grid instead of the sample's own units (ALL 4
# CI_MODELS get this row, not just causal_forest - see
# ci_interval_palette/ci_location_palette above for the full breakdown).
# Split back into an estimator (`model`) plus two *independent* aesthetic
# columns - `interval` (which CI: bootstrap vs. CF's own variance) and
# `location` (evaluated where: the sample's own points vs. the query grid).
# `model` isn't a plotted dimension below - the poster doesn't care about
# differences between the 4 CI_MODELS, so ci_summary pools runs across all
# of them: each point is the bootstrap's (or CF variance estimate's)
# coverage averaged over causal_forest/dr_random_forest/dr_oracle/
# dr_semi_oracle together, not per-estimator. That pooling is a no-op for
# "CF variance estimate", since causal_forest is the only model that ever
# produces that row.
metrics_ci <- metrics_all %>%
  filter(study_name == "confidence_intervals/continuous") %>%
  select(
    scenario,
    n,
    CI_sf,
    run,
    model,
    marginal_coverage,
    simultaneous_coverage
  ) %>%
  mutate(
    interval = if_else(
      model == "causal_forest_inbuilt",
      "CF variance estimate",
      "Half-sample bootstrap"
    ),
    location = if_else(grepl("_grid$", model), "Grid", "Within-sample"),
    model = sub("_inbuilt$|_grid$", "", model)
  ) %>%
  filter(
    scenario %in% c(1, 3, 8, 9) # the poster's usual 4-scenario subset
  ) %>%
  apply_labels(SS_SCENARIO_LABELS) %>%
  mutate(
    n = factor(n, levels = c(500, 1000)),
    CI_sf = factor(CI_sf, levels = seq(0.05, 0.5, 0.05)),
    interval = factor(interval, levels = names(ci_interval_palette)),
    location = factor(location, levels = names(ci_location_palette))
  )

# marginal_coverage is a proportion averaged over many units *within* a run,
# not Bernoulli at the run level, so it gets the general MCSE formula;
# simultaneous_coverage is a genuine 0/1 per run, so it's listed in
# `binomial` for the exact binomial MCSE instead (see summarise_metrics()'s
# roxygen in R/figures.R)
ci_summary <- summarise_metrics(
  metrics_ci,
  c("scenario", "n", "CI_sf", "interval", "location"),
  cols = c(marg_cov = "marginal_coverage", simul_cov = "simultaneous_coverage"),
  binomial = "simul_cov"
)

marg_cov_panel <- point_range_plot(
  ci_summary,
  "marg_cov",
  "Marginal coverage",
  x = "CI_sf",
  colour = "interval",
  shape = "location",
  facet_rows = "n",
  facet_cols = "scenario",
  line = TRUE,
  ci_alpha = 0.5,
  hline = 0.95,
  palette = scale_colour_manual(values = ci_interval_palette),
  shape_palette = scale_shape_manual(values = ci_location_palette)
) +
  labs(x = "Subsampling ratio", colour = "Interval", shape = "Evaluated at") +
  theme(legend.position = "bottom")
ggsave(
  "ci_marginal_coverage.png",
  marg_cov_panel,
  path = figdir,
  height = 20,
  width = 30,
  units = "cm"
)

# doc is undecided whether marginal or simultaneous coverage is the better
# lead metric for the poster - both drafted, pick one later
simul_cov_panel <- point_range_plot(
  ci_summary,
  "simul_cov",
  "Simultaneous coverage",
  x = "CI_sf",
  colour = "interval",
  shape = "location",
  facet_rows = "n",
  facet_cols = "scenario",
  line = TRUE,
  ci_alpha = 0.5,
  hline = 0.95,
  palette = scale_colour_manual(values = ci_interval_palette),
  shape_palette = scale_shape_manual(values = ci_location_palette)
) +
  labs(x = "Subsampling ratio", colour = "Interval", shape = "Evaluated at") +
  theme(legend.position = "bottom")
ggsave(
  "ci_simultaneous_coverage.png",
  simul_cov_panel,
  path = figdir,
  height = 20,
  width = 30,
  units = "cm"
)


# Model evaluation ----------
# This study has its OWN scenario numbering - confirmed only 1/4/6/9 are
# present, not the sample-size studies' 1/3/8/9 - matching
# model_evaluation/me_results.qmd's scenario_labels, not SS_SCENARIO_LABELS.
me_scenario_labels <- c(
  `1` = "Null",
  `4` = "Simple",
  `6` = "Interaction",
  `9` = "Non-linear"
)

# me_metrics.R names every score column <score_type>_<fold_type>_<pipeline>.
# Restricting to the fixed-pi ("05") columns is the doc's own choice ("I only
# want to show the metrics where the propensity was fixed at 0.5"), and the
# main tree (this table) only ever carries fold_type cv/whole - the fuller
# cv_shared/holdout/whole comparison the doc describes lives in
# me_strat_metrics.RDS, not yet collected/synced locally.
me_score_cols <- c(
  "infl05_cv_xgb",
  "infl05_whole_xgb",
  "infl05_cv_automl",
  "infl05_whole_automl",
  "dr05_cv_xgb",
  "dr05_whole_xgb",
  "dr05_cv_automl",
  "dr05_whole_automl"
)

me_long <- metrics_all %>%
  filter(study_name == "model_evaluation") %>%
  select(scenario, n, run, model, true_pehe, all_of(me_score_cols)) %>%
  pivot_longer(
    all_of(me_score_cols),
    names_to = "proxy_id",
    values_to = "score"
  ) %>%
  mutate(
    # simplified version of me_results.qmd's parse_proxy() - no
    # pi_regime/k_groups axes needed since we already filtered to fixed-pi
    # and dropped the calibration score family
    score_family = sub("^([a-z]+)05_.*$", "\\1", proxy_id),
    fold_type = sub("^[a-z]+05_(.*)_(xgb|automl)$", "\\1", proxy_id),
    pipeline = sub("^.*_(xgb|automl)$", "\\1", proxy_id)
  )

# per (scenario, n, run, proxy): which candidate that proxy would have
# picked, and how far its true PEHE is from the best of the 9 - same
# derivation as me_results.qmd's `selection` chunk
me_selection <- me_long %>%
  filter(!is.na(score), !is.na(true_pehe)) %>%
  group_by(scenario, n, run, score_family, fold_type, pipeline) %>%
  filter(n() >= 2) %>%
  summarise(
    pick_rank = rank(true_pehe, ties.method = "min")[which.min(score)],
    regret = true_pehe[which.min(score)] - min(true_pehe),
    .groups = "drop"
  )

me_selection_summary <- me_selection %>%
  group_by(scenario, n, score_family, fold_type, pipeline) %>%
  summarise(
    mean_regret = mean(regret, na.rm = TRUE),
    mcse_regret = sd(regret, na.rm = TRUE) / sqrt(sum(!is.na(regret))),
    .groups = "drop"
  ) %>%
  mutate(
    scenario = factor(
      me_scenario_labels[as.character(scenario)],
      levels = unname(me_scenario_labels)
    ),
    n = factor(n, levels = sort(unique(n))),
    pipeline_lab = recode(pipeline, xgb = "XGBoost", automl = "AutoML"),
    score_lab = recode(score_family, infl = "Influence", dr = "DR risk"),
    proxy_combo = factor(
      paste(score_lab, pipeline_lab, sep = " / "),
      levels = names(proxy_combo_palette)
    ),
    fold_type = factor(
      recode(fold_type, cv = "CV", whole = "Whole"),
      levels = c("CV", "Whole")
    )
  )

# point_range_plot() has no shape aesthetic (the doc wants crossfitting
# regime as shape, on top of the pipeline x score-family colour), so this one
# is a small bespoke plot rather than a further point_range_plot() special case
me_regret_panel <- ggplot(
  me_selection_summary,
  aes(
    x = n,
    y = mean_regret,
    colour = proxy_combo,
    shape = fold_type,
    group = interaction(proxy_combo, fold_type),
    ymin = mean_regret - qnorm(0.975) * mcse_regret,
    ymax = mean_regret + qnorm(0.975) * mcse_regret
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(position = position_dodge(width = 0.5), linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(
    position = position_dodge(width = 0.5),
    linewidth = 0.3,
    alpha = 0.5
  ) +
  facet_grid(cols = vars(scenario), scales = "free_y") +
  scale_colour_manual(values = proxy_combo_palette) +
  drf_theme() +
  theme(legend.position = "bottom") +
  labs(
    y = "Mean regret (excess true PEHE)",
    x = "Sample size",
    colour = "Score / surrogate model",
    shape = "Data splitting"
  )
ggsave(
  "me_regret.png",
  me_regret_panel,
  path = figdir,
  height = 20,
  width = 30,
  units = "cm"
)

# "true rankings, the axis way" - the doc's own brainstorm: does the true-rank
# distribution look different across the pipeline/score-family margins?
# Adapts me_results.qmd's rank_dist/rank_plot chunk (pick_rank, viridis fill -
# kept for this one plot since rank is ordinal, not a fit for a qualitative
# palette), but marginalises pipeline and score_family as two separate facet
# columns instead of one composite label. Aggregated over n and fold_type to
# keep this exploratory version to one plot; a draft, not a final choice.
me_rank_axis <- me_selection %>%
  mutate(
    pipeline_lab = recode(pipeline, xgb = "XGBoost", automl = "AutoML"),
    score_lab = recode(score_family, infl = "Influence", dr = "DR risk"),
    scenario = factor(
      me_scenario_labels[as.character(scenario)],
      levels = unname(me_scenario_labels)
    )
  ) %>%
  pivot_longer(
    c(pipeline_lab, score_lab),
    names_to = "axis",
    values_to = "level"
  ) %>%
  mutate(
    axis = recode(
      axis,
      pipeline_lab = "Nuisance pipeline",
      score_lab = "Score family"
    )
  ) %>%
  count(scenario, axis, level, pick_rank, name = "runs") %>%
  group_by(scenario, axis, level) %>%
  mutate(prop = runs / sum(runs)) %>%
  ungroup()

me_rank_panel <- ggplot(
  me_rank_axis,
  aes(x = level, y = prop, fill = factor(pick_rank))
) +
  geom_col() +
  facet_grid(rows = vars(scenario), cols = vars(axis), scales = "free_x") +
  scale_fill_viridis_d(direction = -1) +
  drf_theme() +
  labs(
    title = "True rankings of selected candidates, marginalised by design axis",
    subtitle = "1 (best) - 9 (worst), fixed-pi DR risk / influence scores only",
    y = "Proportion of runs",
    x = NULL,
    fill = "True rank of pick"
  )
ggsave(
  "me_true_rankings_axis.png",
  me_rank_panel,
  path = figdir,
  height = 20,
  width = 30,
  units = "cm"
)


# Subgroup validation ----------
# validation/continuous/cts_val_metrics.RDS is a named list of 4
# differently-shaped tibbles (see validation/continuous/cts_val_metrics.R),
# so collect_all_metrics.R excludes it from all_metrics.RDS (it only combines
# flat per-run tibbles) - read directly from the study's own result file,
# same as validation/continuous/cts_val_results.qmd.
val_metrics <- readRDS(file.path(
  dirname(path),
  "results",
  "validation",
  "continuous",
  "cts_val_metrics.RDS"
))

# "there's only one model and one subgroup being defined" (planning doc) -
# dr_random_forest only, matching cts_val_results.qmd's precedent
subgroups <- val_metrics$subgroups %>%
  filter(model == "dr_random_forest") %>%
  select(-scenario) %>% # single scenario in this study - not a useful facet
  apply_labels()

# proportion of successful stage 2 W x subgroup interaction tests - already
# built once in cts_val_results.qmd's `subgroup_prop` chunk; ported here
# largely unchanged, just re-styled with the poster's threshold_palette
# instead of drf_scale()/rcartocolor
subgroups_long <- subgroups %>%
  pivot_longer(
    c(top_pval, bottom_pval),
    names_to = "which",
    values_to = "pval"
  ) %>%
  mutate(
    which = recode(which, top_pval = "Top 10%", bottom_pval = "Bottom 10%"),
    which = factor(which, levels = c("Top 10%", "Bottom 10%"))
  )

subgroup_prop <- subgroups_long %>%
  tidyr::crossing(threshold = c(0.05, 0.1)) %>%
  mutate(
    is_sig = pval < threshold,
    threshold = factor(
      threshold,
      levels = c(0.05, 0.1),
      labels = names(threshold_palette)
    )
  ) %>%
  group_by(which, interim_prop, threshold) %>%
  summarise(
    runs = sum(!is.na(is_sig)),
    prop_sig = mean(is_sig, na.rm = TRUE),
    mcse = sqrt(prop_sig * (1 - prop_sig) / runs),
    .groups = "drop"
  )

# NB: only an 8-row local smoke test (2 interim points, 1 scenario) as of
# this draft - runs, but won't show a meaningful curve until the full HPC
# results are synced down
subgroup_panel <- ggplot(
  subgroup_prop,
  aes(
    x = interim_prop,
    y = prop_sig,
    colour = threshold,
    ymin = prop_sig - qnorm(0.975) * mcse,
    ymax = prop_sig + qnorm(0.975) * mcse
  )
) +
  geom_ribbon(aes(fill = threshold), alpha = 0.15, colour = NA) +
  geom_line() +
  geom_point(size = 1.5) +
  facet_grid(cols = vars(which)) +
  scale_x_continuous(breaks = sort(unique(subgroup_prop$interim_prop))) +
  scale_colour_manual(values = threshold_palette) +
  scale_fill_manual(values = threshold_palette) +
  drf_theme(rotate_x = TRUE) +
  labs(
    title = "Proportion of successful stage 2 interaction tests",
    x = "Stage 1 proportion",
    y = "Proportion with p < threshold",
    colour = "Threshold",
    fill = "Threshold"
  )
ggsave(
  "val_subgroup_power.png",
  subgroup_panel,
  path = figdir,
  height = 12,
  width = 20,
  units = "cm"
)

top_subgroup_plot <- subgroup_prop %>%
  filter(which == "Bottom 10%") %>%
  ggplot(aes(
    x = interim_prop,
    y = prop_sig,
    colour = threshold,
    ymin = prop_sig - qnorm(0.975) * mcse,
    ymax = prop_sig + qnorm(0.975) * mcse
  )) +
  geom_ribbon(aes(fill = threshold), alpha = 0.15, colour = NA) +
  geom_line() +
  geom_point(size = 1.5) +
  facet_grid(cols = vars(which)) +
  scale_x_continuous(breaks = sort(unique(subgroup_prop$interim_prop))) +
  scale_colour_manual(values = threshold_palette) +
  scale_fill_manual(values = threshold_palette) +
  drf_theme() +
  theme(legend.position = "bottom") +
  labs(
    x = "Stage 1 proportion",
    y = "Proportion with p < threshold",
    colour = "Threshold",
    fill = "Threshold"
  )

ggsave(
  "val_top_power.png",
  top_subgroup_plot,
  path = figdir,
  height = 20,
  width = 30,
  units = "cm"
)
