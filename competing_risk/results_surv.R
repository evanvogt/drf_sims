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

metrics <- readRDS(paste0(res_path, "new_format/metrics_competing_risks_tidy.RDS"))

metrics <- metrics %>%
  mutate(
    scenario = factor(
      recode(scenario,
             `1` = "No HTE",
             `2` = "HTE on EOI",
             `3` = "HTE on CE",
             `4` = "HTE on both"
      ),
      levels = c("No HTE", "HTE on EOI", "HTE on CE", "HTE on both")
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
  ggplot(aes(x = n, y = bias, color = model, fill = model)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~scenario) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(title = "Average bias in RMST CATEs",
       y = "Bias",
       x = "Sample size")


# MSE
metrics %>%
  filter(mse < 5) %>%
  ggplot(aes(x = n, y = mse, color = model, fill = model)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~scenario) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(title = "MSE in RMST CATEs",
       y = "MSE",
       x = "Sample size")

# Correlation
metrics %>%
  filter(!scenario %in% c("No HTE", "HTE on CE")) %>%
  ggplot(aes(x = n, y = corr, color = model, fill = model)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~scenario) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(title = "Correlation for RMST CATEs",
       y = expression("Pearson's " * rho),
       x = "Sample size")

