##########
# title: plots for the interim-analysis validation study
##########
# Display labels and palette come from R/figures.R instead of the inline
# recode()/scale_*_paletteer_d() calls repeated in every plot below.
#
# interim_prop is a numeric x-axis here, not a facet: there are 11 of them now
# (cts_val_config.R), and the quantity of interest in every summary below is one
# number per interim point, which reads as a curve. cts_val_results.qmd carries
# the distribution-level versions.

library(here)
library(ggplot2)
library(paletteer)
library(dplyr)
library(tidyr)
library(purrr)
source(here("R", "figures.R"))
source(here("validation/continuous/cts_val_config.R"))

metrics <- readRDS(file.path(study$res_path, "cts_val_metrics.RDS"))

MEASURE_LABELS <- c(tevim = "TE-VIM", shap = "TreeSHAP")

tidy_subgroups <- metrics$subgroups %>%
  select(-scenario) %>%
  apply_labels() %>%
  pivot_longer(c(top_pval, bottom_pval), names_to = "which", values_to = "pval") %>%
  mutate(which = recode(which, top_pval = "Top 10%", bottom_pval = "Bottom 10%"),
         which = factor(which, levels = c("Top 10%", "Bottom 10%")))

interim_breaks <- sort(unique(tidy_subgroups$interim_prop))

# Distribution of subgroup p-values, both subgroups. bottom_pval is meaningful
# for the first time here - it used to be the intercept's p-value.
ggplot(tidy_subgroups, aes(x = pval, colour = interim_prop, group = interim_prop)) +
  stat_ecdf(geom = "step", na.rm = TRUE) +
  geom_vline(xintercept = 0.05, linetype = "dashed") +
  facet_grid(model ~ which) +
  scale_colour_viridis_c() +
  drf_theme() +
  labs(
    title = "Cumulative distribution of chunk-2 subgroup p-values",
    x = "p-value (W x subgroup interaction)", y = "Cumulative proportion",
    colour = "Interim\nproportion"
  )

# Proportion of significant (p<0.05) by model, subgroup and interim proportion
tidy_subgroups %>%
  mutate(is_sig = pval < 0.05) %>%
  group_by(model, which, interim_prop) %>%
  summarise(prop_sig = mean(is_sig, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = interim_prop, y = prop_sig, colour = model)) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_line() +
  geom_point(size = 1.5) +
  facet_wrap(~which) +
  scale_x_continuous(breaks = interim_breaks) +
  drf_scale() +
  drf_theme(rotate_x = TRUE) +
  labs(
    title = "Proportion of subgroups validated in later population",
    x = "Interim proportion", y = "Proportion p<0.05"
  )

# Boxplots of CATE variance, early vs. late chunk
tidy_variances <- metrics$variances %>%
  select(-scenario) %>%
  rename(stage1 = vt1, stage2 = vt2) %>%
  apply_labels() %>%
  pivot_longer(cols = c(stage1, stage2), names_to = "stage", values_to = "variance") %>%
  mutate(stage = recode(stage, stage1 = "Chunk 1", stage2 = "Chunk 2"))

ggplot(tidy_variances, aes(x = factor(interim_prop), y = variance, fill = stage)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  facet_wrap(~model, ncol = 1) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  drf_theme(rotate_x = TRUE) +
  labs(
    title = "CATE Variance in early vs. late stages",
    x = "Interim proportion", y = "Variance", fill = "Chunk"
  )

# Variable importance: both measures now (TE-VIM and surrogate TreeSHAP)
tidy_var_imps <- metrics$var_imps %>%
  apply_labels() %>%
  mutate(measure = label_factor(measure, MEASURE_LABELS))

# Mean rank change per covariate
tidy_var_imps %>%
  group_by(variables, model, measure, interim_prop) %>%
  summarise(mean_abs_change = mean(abs(diff), na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = interim_prop, y = mean_abs_change, colour = variables)) +
  geom_line() +
  facet_grid(model ~ measure) +
  scale_x_continuous(breaks = interim_breaks) +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  drf_theme(rotate_x = TRUE) +
  labs(
    title = "Mean Absolute Change in Variable Importance Rank",
    x = "Interim proportion", y = "Mean |Rank change|", colour = "Covariate"
  )

# Average importance rank per variable, by chunk
tidy_var_imps %>%
  group_by(variables, model, measure, interim_prop) %>%
  summarise(
    mean_vi1 = mean(vi1, na.rm = TRUE),
    mean_vi2 = mean(vi2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(mean_vi1, mean_vi2), names_to = "chunk", values_to = "mean_vi") %>%
  mutate(chunk = recode(chunk, mean_vi1 = "Chunk 1", mean_vi2 = "Chunk 2")) %>%
  ggplot(aes(x = interim_prop, y = mean_vi, colour = variables, linetype = chunk)) +
  geom_line() +
  facet_grid(model ~ measure) +
  scale_x_continuous(breaks = interim_breaks) +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  drf_theme(rotate_x = TRUE) +
  labs(
    title = "Average Variable Importance Across Chunks",
    x = "Interim proportion", y = "Mean Importance Rank",
    colour = "Covariate", linetype = "Chunk"
  )

# Does chunk 1's top-ranked covariate still interact with treatment in chunk 2?
metrics$top_var %>%
  apply_labels() %>%
  mutate(measure = label_factor(measure, MEASURE_LABELS)) %>%
  pivot_longer(c(p_cts, p_split), names_to = "form", values_to = "pval") %>%
  mutate(form = recode(form, p_cts = "Continuous W x X", p_split = "Median split")) %>%
  group_by(model, measure, form, interim_prop) %>%
  summarise(prop_sig = mean(pval < 0.05, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = interim_prop, y = prop_sig, colour = model)) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_line() +
  geom_point(size = 1.5) +
  facet_grid(form ~ measure) +
  scale_x_continuous(breaks = interim_breaks) +
  ylim(0, 1) +
  drf_scale() +
  drf_theme(rotate_x = TRUE) +
  labs(
    title = "Top interim covariate: does its treatment interaction hold up?",
    x = "Interim proportion", y = "Proportion p<0.05"
  )
