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
metrics_summary <- metrics %>%
  group_by(scenario, n, CI_sf, model) %>%
  summarise(
    mean_marg_cov = mean(marginal_coverage, na.rm = T),
    mcse_marg_cov = sd(marginal_coverage, na.rm = T)/sqrt(n()),
    mean_simul_cov = mean(simultaneous_coverage, na.rm = T),
    mcse_simul_cov = sd(simultaneous_coverage, na.rm = T)/sqrt(n()),
    mean_ci_len = mean(mean_ci_length, na.rm = T),
    mcse_ci_len = sd(mean_ci_length, na.rm = T)/sqrt(n()),
    .groups = "drop"
  )


metrics_summary %>%
  filter((scenario %in% c(1, 3, 5, 7, 9)) & (n == 1000) & model != "causal_forest_inbuilt") %>%
  ggplot(aes(x = CI_sf, y = mean_marg_cov, colour = model)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_marg_cov - mcse_marg_cov,
                    ymax = mean_marg_cov + mcse_marg_cov)
                ) +
  facet_grid(~scenario) +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal()


metrics_summary %>%
  filter((scenario %in% c(1, 3, 5, 7, 9)) & (n == 1000) & model != "causal_forest_inbuilt") %>%
  ggplot(aes(x = CI_sf, y = mean_simul_cov, colour = model)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_simul_cov - mcse_simul_cov,
                    ymax = mean_simul_cov + mcse_simul_cov)
  ) +
  facet_grid(~scenario) +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal()

metrics_summary %>%
  filter((scenario %in% c(1, 3, 5, 7, 9)) & (n == 1000) & model != "causal_forest_inbuilt") %>%
  ggplot(aes(x = CI_sf, y = mean_ci_len, colour = model)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_ci_len - mcse_ci_len,
                    ymax = mean_ci_len + mcse_ci_len)
  ) +
  facet_grid(~scenario) +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal()
