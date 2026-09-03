##########
# title: figures for the crossfitting CI pilot
##########
rm(list = ls())
# libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(paletteer)
library(here)
library(patchwork)
library(purrr)
library(readr)

# functions - pulls in ci_variant_levels/ci_variant_labels and (transitively,
# via cf_metrics.R) family_labels
source(here("crossfitting", "confidence_intervals", "cf_ci_metrics.R"))

# paths - a separate results tree from crossfitting/'s, see README
path <- here()
res_path <- file.path(dirname(path), "results", "crossfitting_ci")
fig_path <- file.path(dirname(path), "results", "thesis_figures", "crossfitting_ci")
dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)

metrics <- readRDS(file.path(res_path, "cf_ci_metrics.RDS"))

# tidy up - no `set` column here: training-sample tau only, no interval is
# computed for tau_test (see cf_ci_analysis.R)
metrics_all <- metrics %>%
  mutate(
    scenario = factor(
      case_when(scenario == 1 ~ "Null",
                scenario == 6 ~ "Interaction",
                scenario == 9 ~ "Non-linear"),
      levels = c("Null", "Interaction", "Non-linear")),
    family = factor(recode(family, !!!family_labels), levels = unname(family_labels)),
    variant = factor(recode(variant, !!!ci_variant_labels),
                     levels = unname(ci_variant_labels[ci_variant_levels])),
    ci_method = factor(recode(ci_method, !!!ci_method_labels),
                       levels = unname(ci_method_labels[ci_method_levels]))
  ) %>%
  droplevels()

# every figure below the method comparison is about the half-sample bootstrap,
# the one interval available for all 12 arms. The other two methods only exist
# for the 7 whole-sample/OOB arms, so mixing them in would put three rows of a
# different kind next to nine of one kind on the same axis.
metrics <- filter(metrics_all, ci_method == ci_method_labels[["half_boot"]])

# per variant summaries
# mcse denominator is sqrt(sum(!is.na(x))), matching the non-NA count of that
# column - not sqrt(n()), the whole group's row count regardless of NAs.
# interval_metrics() (R/metrics.R) has no internal na.rm, so one failed CI
# bound in a run makes marginal_coverage/mean_ci_length NA for that whole
# run - exactly the case this denominator needs to handle correctly.
metrics_summary <- metrics %>%
  group_by(scenario, family, variant) %>%
  summarise(
    mean_marg_cov = mean(marginal_coverage, na.rm = T),
    mcse_marg_cov = sd(marginal_coverage, na.rm = T) / sqrt(sum(!is.na(marginal_coverage))),
    # simultaneous_coverage is a genuine 0/1 indicator per run (whole-band
    # coverage), so it gets the exact binomial MCSE, not the general one
    # marginal_coverage above uses (that's a proportion over many units
    # within a run, not Bernoulli at the run level)
    mean_simul_cov = mean(simultaneous_coverage, na.rm = T),
    mcse_simul_cov = sqrt(mean_simul_cov * (1 - mean_simul_cov) / sum(!is.na(simultaneous_coverage))),
    mean_ci_len = mean(mean_ci_length, na.rm = T),
    mcse_ci_len = sd(mean_ci_length, na.rm = T) / sqrt(sum(!is.na(mean_ci_length))),
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
    mean_corr = mean(corr, na.rm = T),
    mcse_corr = sd(corr, na.rm = T) / sqrt(sum(!is.na(corr))),
    mean_time_nuisance = mean(time_nuisance, na.rm = T),
    mean_time_stage2 = mean(time_stage2, na.rm = T),
    mean_time_boot = mean(time_boot, na.rm = T),
    mean_time_total = mean(time_nuisance + time_stage2 + time_boot, na.rm = T),
    .groups = "drop"
  )

# helper: the house summary plot, points with a 95% CI (mean +/- qnorm(1 -
# alpha/2) x MCSE) error bar - see continuous/cts_results.R's summary_plot for
# why this isn't a raw +/- 1x MCSE (~68% coverage, not 95%). Unlike
# cf_results.R's version there's no `set` (train/test) dimension to dodge on
# - one sample per arm/run - so colour by variant instead (redundant with x
# under faceting, so the legend is dropped). hline is optional: not every
# metric here has a natural 0/0.95 reference.
summary_plot <- function(df, mean_col, mcse_col, title, ylab, hline = NULL, alpha = 0.05) {
  z <- qnorm(1 - alpha / 2)
  df <- df %>%
    mutate(est = .data[[mean_col]],
           lo = .data[[mean_col]] - z * .data[[mcse_col]],
           hi = .data[[mean_col]] + z * .data[[mcse_col]])
  p <- ggplot(df, aes(x = variant, y = est, colour = variant, ymin = lo, ymax = hi))
  if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, linetype = "dashed")
  p +
    geom_point(size = 2) +
    geom_errorbar(linewidth = 0.3, width = 0.3) +
    facet_wrap(vars(family, scenario), nrow = n_distinct(df$family), scales = "free") +
    scale_colour_paletteer_d("rcartocolor::Safe") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(title = title, y = ylab, x = "Crossfitting procedure")
}

# --- marginal coverage: the headline -----------------------------------
marg_cov_plot <- metrics %>%
  ggplot(aes(x = variant, y = marginal_coverage, colour = variant)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_wrap(vars(family, scenario), nrow = n_distinct(metrics$family), scales = "free") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(title = "Marginal coverage of the half-sample bootstrap interval",
       y = "Marginal coverage", x = "Crossfitting procedure")
ggsave("cf_ci_marg_cov_all.png", path = fig_path, width = 21, height = 15, units = "cm")

marg_cov_sum_plot <- summary_plot(metrics_summary, "mean_marg_cov", "mcse_marg_cov",
                                  "Mean marginal coverage", "Mean marginal coverage",
                                  hline = 0.95)
ggsave("cf_ci_marg_cov_summary.png", plot = marg_cov_sum_plot, path = fig_path,
       width = 21, height = 15, units = "cm")

# --- simultaneous coverage --------------------------------------------------
# per-run 0/1 indicator - a boxplot is uninformative, so only the
# mean-proportion point/errorbar plot is produced
simul_cov_sum_plot <- summary_plot(metrics_summary, "mean_simul_cov", "mcse_simul_cov",
                                   "Mean simultaneous coverage", "Mean simultaneous coverage",
                                   hline = 0.95)
ggsave("cf_ci_simul_cov_summary.png", plot = simul_cov_sum_plot, path = fig_path,
       width = 21, height = 15, units = "cm")

# --- CI width ----------------------------------------------------------
width_plot <- metrics %>%
  ggplot(aes(x = variant, y = mean_ci_length, colour = variant)) +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_wrap(vars(family, scenario), nrow = n_distinct(metrics$family), scales = "free") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(title = "Width of the half-sample bootstrap interval",
       y = "Mean interval length", x = "Crossfitting procedure")
ggsave("cf_ci_width_all.png", path = fig_path, width = 21, height = 15, units = "cm")

width_sum_plot <- summary_plot(metrics_summary, "mean_ci_len", "mcse_ci_len",
                               "Mean interval width", "Mean interval length")
ggsave("cf_ci_width_summary.png", plot = width_sum_plot, path = fig_path,
       width = 21, height = 15, units = "cm")

# --- interval method comparison (OOB arms only) -----------------------------
# The one place all three interval methods exist for the SAME arm, and the
# pilot's sharpest question: does an expensive bootstrap buy anything over the
# variance grf hands back for free?
#
# Read the two panels together, not separately. grf_normal is a POINTWISE
# interval - it is not built to control simultaneous coverage and is expected to
# score near zero on it, so a fair reading of it is the marginal panel. The two
# bootstrap bands are both simultaneous and come off identical refits, differing
# only in whether in-half units are scored; half_boot_out takes its supremum over
# ~n/2 rather than n units, so it is expected to run narrower.
method_df <- metrics_all %>%
  filter(variant %in% ci_variant_labels[c("scf_oob", "scf_oob_t", "oob_oob", "oob_oob_s",
                                          "oob_oob_manual", "cf_full_oob", "cf_default")]) %>%
  droplevels() %>%
  group_by(scenario, family, variant, ci_method) %>%
  summarise(mean_marg_cov = mean(marginal_coverage, na.rm = T),
            mcse_marg_cov = sd(marginal_coverage, na.rm = T) / sqrt(sum(!is.na(marginal_coverage))),
            mean_ci_len = mean(mean_ci_length, na.rm = T),
            mcse_ci_len = sd(mean_ci_length, na.rm = T) / sqrt(sum(!is.na(mean_ci_length))),
            .groups = "drop")

# same 95% CI convention as summary_plot() above - z x MCSE, not raw MCSE
method_z <- qnorm(0.975)

method_cov_panel <- method_df %>%
  ggplot(aes(x = variant, y = mean_marg_cov, colour = ci_method,
             ymin = mean_marg_cov - method_z * mcse_marg_cov,
             ymax = mean_marg_cov + method_z * mcse_marg_cov)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  geom_errorbar(linewidth = 0.3, width = 0.3, position = position_dodge(width = 0.5)) +
  facet_wrap(vars(family, scenario), nrow = n_distinct(method_df$family), scales = "free") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Mean marginal coverage", x = NULL, colour = "Interval method")

method_width_panel <- method_df %>%
  ggplot(aes(x = variant, y = mean_ci_len, colour = ci_method,
             ymin = mean_ci_len - method_z * mcse_ci_len,
             ymax = mean_ci_len + method_z * mcse_ci_len)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  geom_errorbar(linewidth = 0.3, width = 0.3, position = position_dodge(width = 0.5)) +
  facet_wrap(vars(family, scenario), nrow = n_distinct(method_df$family), scales = "free") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Mean interval length", x = "Whole-sample / OOB arm", colour = "Interval method")

method_plot <- (method_cov_panel / method_width_panel) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title = "Three interval methods on the whole-sample / OOB arms",
    subtitle = "grf's variance is pointwise; both bootstrap bands are simultaneous and share their refits")
ggsave("cf_ci_method_comparison.png", plot = method_plot, path = fig_path,
       width = 24, height = 26, units = "cm")

# --- point-estimate accuracy: bias, MSE, correlation ------------------------
# these are per-arm, not per-method, so the half_boot filter above just picks one
# row per arm. Under the same seed these arms are bit-identical to the production
# study's, so they should reproduce cf_results.R's numbers exactly - not this
# report's focus, but the cheapest possible check that the pilot is wired to the
# same estimators
bias_plot <- metrics %>%
  ggplot(aes(x = variant, y = bias, colour = variant)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_wrap(vars(family, scenario), nrow = n_distinct(metrics$family), scales = "free") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(title = "Bias of CATE point estimates by crossfitting procedure",
       y = "Average bias (estimate - truth)", x = "Crossfitting procedure")
ggsave("cf_ci_bias_all.png", path = fig_path, width = 21, height = 15, units = "cm")

mse_plot <- metrics %>%
  ggplot(aes(x = variant, y = mse, colour = variant)) +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_wrap(vars(family, scenario), nrow = n_distinct(metrics$family), scales = "free") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(title = "MSE of CATE point estimates by crossfitting procedure",
       y = "MSE", x = "Crossfitting procedure")
ggsave("cf_ci_mse_all.png", path = fig_path, width = 21, height = 15, units = "cm")

corr_df <- filter(metrics, scenario != "Null")

corr_plot <- corr_df %>%
  ggplot(aes(x = variant, y = corr, colour = variant)) +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  # free_x only: correlation is already on a bounded, comparable -1..1
  # scale, so y stays shared - only the sparse variant axis needs freeing
  facet_wrap(vars(family, scenario), nrow = n_distinct(corr_df$family), scales = "free_x") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(title = "Correlation between estimated and true CATEs",
       y = "Pearson correlation", x = "Crossfitting procedure")
ggsave("cf_ci_corr_all.png", path = fig_path, width = 21, height = 15, units = "cm")

# --- cost ----------------------------------------------------------------
# what nuisances, the final model, and the bootstrap itself each cost
time_df <- metrics_summary %>%
  select(scenario, family, variant, mean_time_nuisance, mean_time_stage2, mean_time_boot) %>%
  pivot_longer(c(mean_time_nuisance, mean_time_stage2, mean_time_boot),
               names_to = "stage", values_to = "seconds") %>%
  mutate(stage = recode(stage,
                        mean_time_nuisance = "Stage 1 (nuisances)",
                        mean_time_stage2 = "Stage 2 (final model)",
                        mean_time_boot = "Bootstrap CI"))

time_plot <- time_df %>%
  ggplot(aes(x = variant, y = seconds, fill = stage)) +
  geom_col() +
  facet_wrap(vars(family, scenario), nrow = n_distinct(time_df$family), scales = "free") +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Mean runtime per replicate",
       y = "Seconds", x = "Crossfitting procedure", fill = "Stage")
ggsave("cf_ci_runtime.png", path = fig_path, width = 21, height = 15, units = "cm")

# --- coverage against cost ---------------------------------------------
# the decision plot: does the bootstrap buy closer-to-nominal coverage, and
# at what cost, across arms?
coverage_cost_plot <- metrics_summary %>%
  ggplot(aes(x = mean_time_total, y = mean_marg_cov, colour = variant, shape = family)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_point(size = 3) +
  facet_wrap(~scenario, scales = "free") +
  scale_x_log10() +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(title = "Marginal coverage against runtime",
       subtitle = "Closer to the dashed line is better",
       y = "Mean marginal coverage", x = "Mean runtime per replicate (s, log scale)",
       colour = "Procedure", shape = "Estimator")
ggsave("cf_ci_coverage_cost.png", path = fig_path, width = 21, height = 15, units = "cm")

# --- headline table -------------------------------------------------------
headline <- metrics_summary %>%
  select(scenario, family, variant, mean_marg_cov, mcse_marg_cov, mean_simul_cov,
         mean_ci_len, mean_bias, mean_ate_bias, mean_rel_ate_bias,
         mean_rel_bias_cate, mean_mse, mean_time_total) %>%
  arrange(family, scenario, desc(mean_marg_cov))

print(headline, n = Inf)
saveRDS(headline, file.path(res_path, "cf_ci_headline.RDS"))
write_csv(headline, file.path(res_path, "cf_ci_headline.csv"))

print("figures written!")
