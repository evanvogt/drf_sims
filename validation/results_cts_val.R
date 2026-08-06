##########
# title: plots for the interim-analysis validation study
##########
# Display labels and palette come from R/figures.R instead of the inline
# recode()/scale_*_paletteer_d() calls repeated in every plot below.

library(here)
library(ggplot2)
library(paletteer)
library(dplyr)
library(tidyr)
library(purrr)
source(here("R", "figures.R"))
source(here("validation/cts_val_config.R"))

metrics <- readRDS(file.path(study$res_path, "cts_val_metrics.RDS"))

tidy_subgroups <- metrics$subgroups %>%
  select(-scenario) %>%
  apply_labels()

# Histogram of top p-values
ggplot(tidy_subgroups, aes(x = top_pval, fill = model)) +
  geom_histogram(bins = 30, alpha = 0.5, position = "identity") +
  facet_wrap(~interim_prop) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  drf_scale() +
  drf_theme() +
  labs(
    title = "Histogram of Top 10% Subgroup p-values",
    x = "p-value (Top 10%)", y = "Count"
  )

# Proportion of significant (p<0.05) by model and interim proportion
tidy_subgroups %>%
  mutate(is_sig = top_pval < 0.05) %>%
  group_by(model, interim_prop) %>%
  summarise(prop_sig = mean(is_sig, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = interim_prop, y = prop_sig, fill = model)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  drf_scale() +
  drf_theme() +
  labs(
    title = "Proportion of subgroups validated in later population",
    x = "Interim proportion", y = "Proportion p<0.05"
  )

# Boxplots of CATE variance, early vs. late chunk
tidy_variances <- metrics$variances %>%
  select(-scenario) %>%
  rename(stage1 = vt1, stage2 = vt2) %>%
  apply_labels() %>%
  pivot_longer(cols = c(stage1, stage2), names_to = "stage", values_to = "variance")

ggplot(tidy_variances, aes(x = stage, y = variance, fill = model, color = model)) +
  geom_boxplot(alpha = 0.7) +
  facet_grid(model ~ interim_prop) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  drf_scale() +
  drf_theme() +
  labs(
    title = "CATE Variance in early vs. late stages",
    x = "", y = "Variance"
  ) +
  theme(legend.position = "none")

# Mean rank change barplot
top_vars <- metrics$var_imps %>%
  apply_labels() %>%
  group_by(variables, model, interim_prop) %>%
  summarise(mean_abs_change = mean(abs(diff), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abs_change))

ggplot(top_vars, aes(x = reorder(variables, mean_abs_change), y = mean_abs_change, fill = model)) +
  geom_col(position = "dodge") +
  coord_flip() +
  facet_wrap(~interim_prop) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  drf_scale() +
  drf_theme() +
  labs(
    title = "Mean Absolute Change in Variable Importance Rank",
    x = "Variable", y = "Mean |Rank change|"
  )

# Lineplot of average importance per variable, by chunk
mean_vi <- metrics$var_imps %>%
  apply_labels() %>%
  group_by(variables, model) %>%
  summarise(
    mean_vi1 = mean(vi1, na.rm = TRUE),
    mean_vi2 = mean(vi2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(mean_vi1, mean_vi2), names_to = "chunk", values_to = "mean_vi")

ggplot(mean_vi, aes(x = variables, y = mean_vi, color = chunk, group = chunk)) +
  geom_point() + geom_line() +
  facet_wrap(~model) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  drf_scale() +
  drf_theme() +
  coord_flip() +
  labs(
    title = "Average Variable Importance Across Chunks",
    x = "Variable", y = "Mean Importance Rank"
  )
