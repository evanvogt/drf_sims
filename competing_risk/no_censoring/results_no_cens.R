#####################
# title: plots for continuous results
####################


# libraries and data
library(ggplot2)
library(paletteer)
library(gridExtra)
library(dplyr)
library(tidyr)
library(purrr)

res_path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/results/"

metrics <- readRDS(paste0(res_path, "new_format/metrics_competing_risks_no_cens_tidy.RDS"))

metrics <- metrics %>%
  mutate(
    scenario = factor(
      recode(scenario,
             `1` = "ATE event 1",
             `2` = "ATE event 2",
             `3` = "HTE event 1",
             `4` = "HTE event 1, ATE event 2",
             `5` = "HTE event 2",
             `6` = "HTE event 2, ATE event 1",
             `7` = "HTE event 1 & 2"),
      levels = c("ATE event 1", "ATE event 2", "HTE event 1", "HTE event 1, ATE event 2", "HTE event 2", "HTE event 2, ATE event 1", "HTE event 1 & 2")
    ),
    n = factor(
      recode(n,
             `250` = "250",
             `500` = "500",
             `1000` = "1000"),
      levels = c("250", "500", "1000")
    ),
    model = factor(
      recode(model,
             causal_survival_forest_sub_dist = "CSF sub-distrubtion",
             causal_survival_forest_cause_spec = "CSF cause-specific"),
      levels = c("CSF sub-distrubtion", "CSF cause-specific"))
  )
# bias
metrics %>%
  filter(abs(bias) < 10) %>%
  ggplot(aes(x = n, y = bias, color = model, fill = model)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~scenario) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "Bias",
       x = "Sample size")


# MSE
metrics %>%
  filter(mse < 20) %>%
  ggplot(aes(x = n, y = mse, color = model, fill = model)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~scenario) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "MSE",
       x = "Sample size")

# Correlation
metrics %>%
  ggplot(aes(x = n, y = corr, color = model, fill = model)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~scenario) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = expression("Pearson's " * rho),
       x = "Sample size")

