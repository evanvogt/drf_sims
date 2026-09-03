##########
# title: results of confidence interval stuff
##########


# libraries
library(ggplot2)
library(paletteer)
library(dplyr)
library(here)

# paths
path <- here()
res_path <- file.path(dirname(here()), "results", "confidence_intervals", "continuous")

# data
metrics <- readRDS(file.path(res_path, "ci_cts_metrics.RDS"))

# per scenario summaries
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

metrics_summary %>%
  filter((scenario %in% c(1, 3, 5, 7, 9)) & (n == 1000) & model != "causal_forest_inbuilt") %>%
  ggplot(aes(x = CI_sf, y = mean_marg_cov, colour = model)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_marg_cov - z * mcse_marg_cov,
                    ymax = mean_marg_cov + z * mcse_marg_cov)
                ) +
  facet_grid(~scenario) +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal()


metrics_summary %>%
  filter((scenario %in% c(1, 3, 5, 7, 9)) & (n == 1000) & model != "causal_forest_inbuilt") %>%
  ggplot(aes(x = CI_sf, y = mean_simul_cov, colour = model)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_simul_cov - z * mcse_simul_cov,
                    ymax = mean_simul_cov + z * mcse_simul_cov)
  ) +
  facet_grid(~scenario) +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal()

metrics_summary %>%
  filter((scenario %in% c(1, 3, 5, 7, 9)) & (n == 1000) & model != "causal_forest_inbuilt") %>%
  ggplot(aes(x = CI_sf, y = mean_ci_len, colour = model)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_ci_len - z * mcse_ci_len,
                    ymax = mean_ci_len + z * mcse_ci_len)
  ) +
  facet_grid(~scenario) +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal()
