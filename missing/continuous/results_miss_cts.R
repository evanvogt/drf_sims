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

metrics <- readRDS(paste0(res_path, "new_format/metrics_continuous_miss_tidy.RDS"))

metrics <- metrics %>%
  select(-c(scenario, n)) %>%
  mutate(
    model = factor(
      recode(model,
             causal_forest = "Causal forest",
             dr_random_forest = "DR-RandomForest",
             dr_superlearner = "DR-SuperLearner"),
      levels = c("Causal forest", "DR-RandomForest", "DR-SuperLearner")
    ),
    miss_type = factor(miss_type,
      levels = c("prognostic", "predictive", "both")
    ),
    miss_prop = factor(
      recode(miss_prop,
             `0.2` = "0.2",
             `0.4` = "0.4"),
      levels = c("0.2", "0.4")
    ),
    miss_method = factor(
      recode(miss_method,
             complete_cases = "complete cases",
             mean_imputation = "mean imputation",
             multiple_imputation = "mulitple imputation"),
      levels = c("complete cases", "mean imputation", "mulitple imputation")
    )
  )
# bias
metrics %>%
  ggplot(aes(x = miss_prop, y = bias, color = model, fill = model)) +
  geom_boxplot(alpha = 0.7) +
  facet_grid(miss_method~miss_type) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(title = "Average bias in continuous CATEs with missing data",
       y = "Bias",
       x = "missingness proportions")


# MSE
metrics %>%
#  filter(mse < 5) %>%
  ggplot(aes(x = miss_prop, y = mse, color = model, fill = model)) +
  geom_boxplot(alpha = 0.7) +
  facet_grid(miss_method~miss_type) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(title = "MSE in continuous CATEs with missing data",
       y = "MSE",
       x = "Missingness proportion")

# Correlation
metrics %>%
  ggplot(aes(x = miss_prop, y = corr, color = model, fill = model)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(miss_method~miss_type) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(title = "Correlation with true CATEs under missingness",
       y = expression("Pearson's " * rho),
       x = "Sample size")

# correct tests 
metrics %>%
  group_by(miss_prop, miss_type, miss_method, model) %>%
  summarise(across(c(BLP_correct), mean), .groups = "drop") %>%
  pivot_longer(
    cols = c(BLP_correct),
    names_to = "Metric",
    values_to = "Proportion"
  ) %>%
  ggplot(aes(x = miss_prop, y = Proportion, color = model, fill = model)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  facet_grid(miss_method~miss_type) +
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
