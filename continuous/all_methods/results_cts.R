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

metrics <- readRDS(paste0(res_path, "new_format/metrics_continuous_tidy.RDS"))

metrics <- metrics %>%
  mutate(
    scenario = factor(
      recode(scenario,
             `1` = "No HTE",
             `3` = "Simple HTE",
             `8` = "Complex HTE",
             .default = "Non-linear HTE"
      ),
      levels = c("No HTE", "Simple HTE", "Complex HTE", "Non-linear HTE")
    ),
    n = factor(
      recode(n,
             `100` = "100",
             `250` = "250",
             `500` = "500",
             `1000` = "1000"),
      levels = c("100", "250", "500", "1000")
    ),
    model = factor(
      recode(model,
             causal_forest = "Causal forest",
             dr_random_forest = "DR-RandomForest",
             dr_superlearner = "DR-SuperLearner"),
      levels = c("Causal forest", "DR-RandomForest", "DR-SuperLearner"))
  )
# bias
metrics %>%
  ggplot(aes(x = n, y = bias, color = model, fill = model)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~scenario) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(title = "Average bias in continuous CATEs",
       y = "Bias",
       x = "Sample size")


# MSE
metrics %>%
  filter(mse < 3) %>%
  ggplot(aes(x = n, y = mse, color = model, fill = model)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~scenario) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(title = "MSE in continuous CATEs",
       y = "MSE",
       x = "Sample size")

# Correlation
metrics %>%
  filter(scenario != "No HTE") %>%
  ggplot(aes(x = n, y = corr, color = model, fill = model)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~scenario) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(title = "MSE in continuous CATEs",
       y = expression("Pearson's " * rho),
       x = "Sample size")

# correct tests 
metrics %>%
  group_by(scenario, n, model) %>%
  summarise(across(c(BLP_correct, indep_correct), mean), .groups = "drop") %>%
  pivot_longer(
    cols = c(BLP_correct, indep_correct),
    names_to = "Metric",
    values_to = "Proportion"
  ) %>%
  ggplot(aes(x = n, y = Proportion, color = model, fill = model)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  facet_grid(rows = vars(Metric), cols = vars(scenario)) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(title = "BLP tests",
       y = "proportion of simulations",
       x = "Sample size")

summary_table <- metrics %>%
  group_by(scenario, n, model) %>%
  summarise(
    mean_BLP_correct = mean(BLP_correct, na.rm = TRUE),
    mean_indep_correct = mean(indep_correct, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_table)
