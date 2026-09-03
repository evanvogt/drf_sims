##########
# title: figures for the thesis chapter - cts miss CI
##########

# libraries
# Labels, palette and figure sizing come from R/figures.R.
library(here)
library(patchwork)
library(purrr)
library(ggridges)
library(scales)
source(here("R", "figures.R"))
# paths
path <- here()
res_path <- file.path(dirname(path), "results", "missing", "ci_example")
fig_path <- file.path(dirname(path), "results", "thesis_figures", "miss_cts_ci")
dir.create(fig_path, showWarnings = F, recursive = T)

metrics <- readRDS(file.path(res_path, "cts_miss_ci_metrics.RDS"))

# tidy up. MISS_SCENARIO_LABELS and the pooling-strategy labels now come from
# R/figures.R; `strategy` is what cts_miss_ci_metrics.R actually writes (this
# script used to name it CI_method, a column the metrics file never had).
metrics <- metrics %>%
  filter(scenario %in% c(1, 2, 4, 5)) %>%
  apply_labels(MISS_SCENARIO_LABELS) %>%
  mutate(n = factor(n, levels = c(100, 250, 500, 1000)))


# per scenario summaries
#
# mcse denominator is sqrt(sum(!is.na(x))), matching the non-NA count of that
# column - not sqrt(n()), the whole group's row count regardless of NAs.
# mcse_mar_cov/mcse_sim_cov/mcse_ci_len previously had NO denominator at all
# (just sd(x)) - those were raw standard deviations plotted as if they were
# standard errors, inflating the error bars by sqrt(n_runs).
metrics_summary <- metrics %>%
  group_by(scenario, n, type, prop, mechanism, method, model, strategy) %>%
  summarise(
    mean_bias = mean(bias, na.rm = T),
    mcse_bias = sd(bias, na.rm = T)/sqrt(sum(!is.na(bias))),
    mean_mse = mean(mse, na.rm = T),
    mcse_mse = sd(mse, na.rm = T)/sqrt(sum(!is.na(mse))),
    mean_mar_cov = mean(marginal_coverage, na.rm = T),
    mcse_mar_cov = sd(marginal_coverage, na.rm = T)/sqrt(sum(!is.na(marginal_coverage))),
    # simultaneous_coverage is a genuine 0/1 indicator per run (whole-band
    # coverage), so it gets the exact binomial MCSE, not the general one
    # marginal_coverage above uses (that's a proportion over many units
    # within a run, not Bernoulli at the run level)
    mean_sim_cov = mean(simultaneous_coverage, na.rm = T),
    mcse_sim_cov = sqrt(mean_sim_cov * (1 - mean_sim_cov)/sum(!is.na(simultaneous_coverage))),
    mean_ci_len = mean(mean_ci_length, na.rm = T),
    mcse_ci_len = sd(mean_ci_length, na.rm = T)/sqrt(sum(!is.na(mean_ci_length))),
    .groups = "drop"
  )

# error bars below are a 95% CI (mean +/- qnorm(0.975) x MCSE), not a raw
# +/- 1x MCSE (~68% coverage) - see continuous/cts_results.R's summary_plot
z <- qnorm(0.975)

# combining everything into a single plot?


# marginal coverage
mar_cov_plot <- metrics %>%
  ggplot(aes(x=strategy, y = marginal_coverage, colour = model)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_wrap(~scenario) +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "Marginal coverage",
       x = "Method")

# summary
mar_cov_sum_plot <- metrics_summary %>%
  ggplot(aes(x = strategy, y = mean_mar_cov, colour = model, ymin = mean_mar_cov - z * mcse_mar_cov, ymax = mean_mar_cov + z * mcse_mar_cov)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(position = position_dodge(width = 0.5), linewidth = 0.3) +
  facet_wrap(~scenario) +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "Marginal Coverage",
       x = "Method")

# simultaneous coverage
sim_cov_plot <- metrics %>%
  ggplot(aes(x=strategy, y = simultaneous_coverage, colour = model)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_wrap(~scenario) +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "Simultaneous coverage",
       x = "Method") +
  theme(axis.text.x = element_text(angle = 90))

# summary
sim_cov_sum_plot <- metrics_summary %>%
  ggplot(aes(x = strategy, y = mean_sim_cov, colour = model, ymin = mean_sim_cov - z * mcse_sim_cov, ymax = mean_sim_cov + z * mcse_sim_cov)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(position = position_dodge(width = 0.5), linewidth = 0.3) +
  facet_wrap(~scenario) +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "Simultaneous coverage",
       x = "Method")

# CI length
ci_len_plot <- metrics %>%
  ggplot(aes(x=strategy, y = mean_ci_length, colour = model)) +
  geom_boxplot(fill = "transparent", outlier.shape = NA) +
  facet_wrap(~scenario) +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "CI length",
       x = "Method")

# summary
ci_len_sum_plot <- metrics_summary %>%
  ggplot(aes(x = strategy, y = mean_ci_len, colour = model, ymin = mean_ci_len - z * mcse_ci_len, ymax = mean_ci_len + z * mcse_ci_len)) +
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
