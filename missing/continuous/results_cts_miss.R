##########
# title: visualising cts missingness results
##########


# libraries
library(ggplot2)
library(paletteer)
library(dplyr)

# paths
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live"

metrics <- readRDS(file.path(path, "results/new_format/metrics_cts_miss_df.RDS"))


# means and mcses df
metrics_summary <- metrics %>%
  group_by(scenario, type, mechanism, method, model) %>%
  summarise(
    mean_bias = mean(bias, na.rm = T),
    mcse_bias = sd(bias, na.rm = T) / sqrt(sum(!is.na(bias))),
    mean_mse = mean(mse, na.rm = T),
    mcse_mse = sd(mse, na.rm = T) / sqrt(sum(!is.na(mse))),
    .groups = "drop"
  )

# bias
metrics %>%
  filter(type == "both") %>%
  filter(scenario == "scenario_5") %>%
  #filter(model %in% c("Causal forest", "DR-RandomForest")) %>%
  ggplot(aes(x = model, y = bias, color = model, fill = model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(alpha = 0.7) +
  facet_grid(mechanism~method) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(colour = "black")) +
  labs(title = "Average bias in continuous CATEs with missing data",
       y = "Bias",
       x = "missing data handling method") +
  theme(axis.text.x = element_blank())

# error bars are a 95% CI (mean +/- qnorm(0.975) x MCSE), not a raw +/- 1x
# MCSE (~68% coverage)
z <- qnorm(0.975)

metrics_summary %>%
  ggplot(aes(x = model, y = mean_bias, color = model, fill = model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_bias - z * mcse_bias,
                    ymax = mean_bias + z * mcse_bias),
                width = 0.2) +
  facet_grid(mechanism ~ method) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(colour = "black")) +
  labs(title = "Mean bias (95% MC CI) in continuous CATEs with missing data",
       y = "Bias",
       x = "missing data handling method") +
  theme(axis.text.x = element_blank())


metrics %>%
  filter(type == "both" & mse < 2) %>%
  ggplot(aes(x = model, y = mse, color = model, fill = model)) +
  geom_boxplot(alpha = 0.7) +
  facet_grid(method~mechanism) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(colour = "black")) +
  labs(title = "Average MSE in continuous CATEs with missing data",
       y = "MSE",
       x = "missingness proportions")
