##########
# title: figures for the thesis chapter - cts outcome confidence intervals
##########

# Labels, palette and figure sizing come from R/figures.R.
library(here)
library(patchwork)
library(purrr)
library(ggridges)
library(scales)
source(here("R", "figures.R"))

# paths
path <- here()
res_path <- file.path(dirname(path), "results", "confidence_intervals", "continuous")
fig_path <- file.path(dirname(path), "results", "thesis_figures", "cts_ci")
dir.create(fig_path, showWarnings = F, recursive = T)

metrics <- readRDS(file.path(res_path, "ci_cts_metrics.RDS"))

metrics <- metrics %>%
  filter(scenario %in% c(1, 3, 8, 9)) %>%
  apply_labels(SS_SCENARIO_LABELS) %>%
  # model labels (including causal_forest_inbuilt) come from apply_labels()
  mutate(
    n = factor(n, levels = c(500, 1000)),
    CI_sf = factor(CI_sf, levels = seq(0.05, 0.5, 0.05))
  )

# mcse denominator is sqrt(sum(!is.na(x))), matching the non-NA count of that
# column - not sqrt(n()), the whole group's row count regardless of NAs.
# interval_metrics() (R/metrics.R) has no internal na.rm, so one failed CI
# bound in a run makes marginal_coverage/mean_ci_length NA for that whole run.
metrics_summary <- metrics %>%
  group_by(scenario, n, CI_sf, model) %>%
  summarise(
    mean_marg_cov = mean(marginal_coverage, na.rm = T),
    mcse_marg_cov = sd(marginal_coverage, na.rm = T)/sqrt(sum(!is.na(marginal_coverage))),
    # simultaneous_coverage is a genuine 0/1 indicator per run (whole-band
    # coverage), so it gets the exact binomial MCSE, not the general one
    # marginal_coverage above uses (that's a proportion over many units
    # within a run, not Bernoulli at the run level)
    mean_simul_cov = mean(simultaneous_coverage, na.rm = T),
    mcse_simul_cov = sqrt(mean_simul_cov * (1 - mean_simul_cov)/sum(!is.na(simultaneous_coverage))),
    mean_ci_len = mean(mean_ci_length, na.rm = T),
    mcse_ci_len = sd(mean_ci_length, na.rm = T)/sqrt(sum(!is.na(mean_ci_length))),
    .groups = "drop"
  )

# error bars are a 95% CI (mean +/- qnorm(0.975) x MCSE), not a raw +/- 1x
# MCSE (~68% coverage) - see continuous/cts_results.R's summary_plot
z <- qnorm(0.975)

mc_sum_plot <- metrics_summary %>%
  ggplot(aes(x = CI_sf, y = mean_marg_cov, colour = model, ymin = mean_marg_cov - z * mcse_marg_cov, ymax = mean_marg_cov + z * mcse_marg_cov)) +
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
  ggplot(aes(x = CI_sf, y = mean_simul_cov, colour = model, ymin = mean_simul_cov - z * mcse_simul_cov, ymax = mean_simul_cov + z * mcse_simul_cov)) +
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
  ggplot(aes(x = CI_sf, y = mean_ci_len, colour = model, ymin = mean_ci_len - z * mcse_ci_len, ymax = mean_ci_len + z * mcse_ci_len)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(position = position_dodge(width = 0.5)) +
  facet_grid(cols = vars(scenario), rows = vars(n), scales = "free_y", axes = "all_x") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "Confidence interval length",
       x = "Subsampling ratio")
ggsave("cts_ci_ci_len_sum.png", path = fig_path, width = 21, height = 15, units = "cm")
