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
    mcse_bias = sd(bias, na.rm = T)/sqrt(n()),
    mean_mse = mean(mse, na.rm = T),
    mcse_mse = sd(mse, na.rm = T)/sqrt(n()),
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
  theme_minimal() +
  labs(title = "Average bias in continuous CATEs with missing data",
       y = "Bias",
       x = "missing data handling method") +
  theme(axis.text.x = element_blank())

metrics_summary %>%
  ggplot(aes(x = model, y = mean_bias, color = model, fill = model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_bias - mcse_bias,
                    ymax = mean_bias + mcse_bias),
                width = 0.2) +
  facet_grid(mechanism ~ method) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(title = "Mean bias (± MCSE) in continuous CATEs with missing data",
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
  theme_minimal() +
  labs(title = "Average MSE in continuous CATEs with missing data",
       y = "MSE",
       x = "missingness proportions")
