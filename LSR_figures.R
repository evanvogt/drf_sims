##############
# LSR metrics
##############


# libraries and data
library(ggplot2)
library(paletteer)
library(gridExtra)
library(dplyr)
library(tidyr)
library(purrr)

res_path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/results/"

metrics_cts <- readRDS(paste0(res_path, "new_format/metrics_continuous_tidy.RDS"))

metrics_bin <- readRDS(paste0(res_path, "new_format/metrics_binary_tidy.RDS"))

res_bin <- metrics_bin %>%
  filter(scenario %in% c(1,3,8,9)) %>%
  mutate(
    scenario = factor(
      recode(scenario,
             `1` = "Null",
             `3` = "Simple",
             `8` = "Complex",
             `9` = "Non-linear"),
      levels = c("Null", "Simple", "Complex", "Non-linear")
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
             dr_superlearner = "DR-SuperLearner",
             dr_oracle = "DR-oracle",
             dr_semi_oracle = "DR-semi-oracle"),
      levels = c("Causal forest", "DR-RandomForest", "DR-SuperLearner", "DR-oracle", "DR-semi-oracle"))
  )
# bias
res_bin %>%
  filter(abs(bias) < 1.4) %>%
  ggplot(aes(x = n, y = bias, color = model, fill = model)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~scenario, nrow = 2) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "Bias",
       x = "Sample size")


# MSE
res_bin %>%
  filter(mse < 0.5) %>%
  ggplot(aes(x = n, y = mse, color = model, fill = model)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~scenario, nrow = 2) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "MSE",
       x = "Sample size")

# removing the 100 scenario to see performance betweeen models better:
res_cts %>%
  filter(mse < 2) %>%
  filter(n != "100") %>%
  ggplot(aes(x = n, y = mse, color = model, fill = model)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~scenario, nrow = 2) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(title = "MSE in continuous CATEs",
       y = "MSE",
       x = "Sample size")

# Correlation
res_bin %>%
  filter(scenario != "Null") %>%
  ggplot(aes(x = n, y = corr, color = model, fill = model)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~scenario) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = expression("Pearson's " * rho),
       x = "Sample size")

# correct tests 
res_bin %>%
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
  labs(title = "proportion of correct BLP tests",
       y = "proportion of simulations",
       x = "Sample size")

# distribution of p-values
res_bin %>%
  pivot_longer(
    cols = c(BLP_p, indep_p),
    names_to = "Test",
    values_to = "p_value"
  ) %>%
  ggplot(aes(x = n, y = p_value, color = model, fill = model)) +
  geom_violin(alpha = 0.7) +
  facet_grid(rows = vars(Test), cols = vars(scenario)) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  coord_flip() +
  labs(title = "P-value distribution from HTE tests",
       y = "p-value",
       x = "Sample size")

# separate plots by test
library(ggridges)
res_bin %>%
  ggplot(aes(x = BLP_p, y = model, color = model, fill = model, group = interaction(n, scenario, model))) +
  geom_density_ridges(alpha = 0.5) +
  facet_grid(rows = vars(n), cols = vars(scenario)) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  xlim(0,1) +
  labs(title = "P-value distribution from BLP test",
       x = "p-value",
       y = "Sample size")

res_cts %>%
  ggplot(aes(x = indep_p, y = model, color = model, fill = model, group = interaction(n, scenario, model))) +
  geom_density_ridges2(alpha = 0.5) +
  facet_grid(rows = vars(n), cols = vars(scenario)) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  xlim(0,1) +
  labs(title = "P-value distribution from permutation test",
       x = "p-value",
       y = "Sample size") +
  theme(axis.text.x = element_text(size = 5))

# not splitting by the model

res_bin %>%
  pivot_longer(
    cols = c(BLP_p, indep_p),
    names_to = "Test",
    values_to = "p_value"
  ) %>%
  mutate(Test = factor(recode(Test, `BLP_p` = "BLP", `indep_p` = "Permutation"), levels = c("BLP", "Permutation"))) %>%
  ggplot(aes(x = p_value, y = n, color = n, fill = n, group = interaction(n, scenario, Test))) +
  geom_density_ridges2(alpha = 0.5) +
  facet_grid(Test~scenario) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  xlim(0,1) +
  labs(x = "p-value",
       y = "Sample size") +
  theme(axis.text.x = element_text(size = 5))

# table of correctness
table_cts <- res_cts %>%
  group_by(scenario, n) %>%
  summarise(
    mean_BLP_correct = mean(BLP_correct, na.rm = TRUE),
    mean_indep_correct = mean(indep_correct, na.rm = TRUE),
    .groups = "drop"
  )

print(table_cts)
print(table_bin)
