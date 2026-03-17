##########
# title: figures for the thesis chapter - bin outcomes
##########

# libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(paletteer)
library(here)
library(patchwork)
library(purrr)
library(ggridges)
library(scales)

# paths
path <- here()
res_path <- file.path(dirname(path), "results", "missing", "binary")
fig_path <- file.path(dirname(path), "results", "thesis_figures", "miss_bin")
dir.create(fig_path, showWarnings = F, recursive = T)

metrics <- readRDS(file.path(res_path, "bin_miss_metrics.RDS"))

metrics <- metrics %>%
  filter((scenario %in% c(1, 2, 4, 5) & type == "both") | (scenario == 1)) %>%
  mutate(scenario = factor(
    case_when(scenario == 1 ~ "Null",
              scenario == 2 ~ "Simple",
              scenario == 4 ~ "Complex",
              scenario == 5 ~ "Non-linear"), levels = c("Null", "Simple", "Complex", "Non-linear")),,
    n = factor(n, levels = c(500)),
    prop = factor(prop, levels = c(0.3)),
    mechanism = factor(
      recode(mechanism,
             MAR = "MAR",
             AUX = "MNAR",
             `AUX-Y` = "MNAR-Y"),
      levels = c("MAR", "MNAR", "MNAR-Y")),
    method = factor(
      recode(method,
             complete_cases = "Complete case",
             mean_imputation = "Single mean imputation",
             missforest = "Single forest-based imputation",
             regression = "Single model-based imputation",
             missing_indicator = "Missing indicators",
             IPW = "IPW",
             multiple_imputation = "Multiple imputation (rf)",
             none = "Inbuilt missingness handling"),
      levels = c("Complete case", "Single mean imputation", "Single forest-based imputation", "Single model-based imputation", "Missing indicators", "IPW", "Multiple imputation (rf)", "Inbuilt missingness handling")
    ),
    model = factor(
      recode(model,
             causal_forest = "Causal forest",
             dr_random_forest = "DR-RandomForest",
             dr_oracle = "DR-oracle",
             dr_semi_oracle = "DR-semi-oracle",
             dr_superlearner = "DR-SuperLearner"),
      levels = c("Causal forest", "DR-RandomForest", "DR-oracle", "DR-semi-oracle", "DR-SuperLearner")
    ))

metrics_summary <- metrics %>%
  group_by(scenario, n, type, prop, mechanism, method, model) %>%
  summarise(
    mean_bias = mean(bias, na.rm = T),
    mcse_bias = sd(bias, na.rm = T)/sqrt(n()),
    mean_mse = mean(mse, na.rm = T),
    mcse_mse = sd(mse, na.rm = T)/sqrt(n()),
    mean_corr = mean(corr, na.rm = T),
    mcse_corr = sd(corr, na.rm = T)/sqrt(n()),
    .groups = "drop"
  )


# bias
bias_sum_plot <-  metrics_summary %>%
  ggplot(aes(x = model, y = mean_bias, colour = method, ymin = mean_bias - mcse_bias, ymax = mean_bias + mcse_bias)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(position = position_dodge(width = 0.5), linewidth = 0.3) +
  facet_grid(rows = vars(mechanism), cols = vars(scenario), scale = "free") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "Bias",
       x = "Model") +
  theme(axis.text.x = element_text(angle = 90))
ggsave("miss_bin_bias.png", path = fig_path, width = 21, height = 15, units = "cm")

mse_sum_plot <- metrics_summary %>%
  ggplot(aes(x = model, y = mean_mse, colour = method, ymin = mean_mse - mcse_mse, ymax = mean_mse + mcse_mse)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(position = position_dodge(width = 0.5), linewidth = 0.3) +
  facet_grid(rows = vars(mechanism), cols = vars(scenario), scales = "free_y") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "MSE",
       x = "Model") +
  theme(axis.text.x = element_text(angle = 90))
ggsave("miss_bin_mse.png", path = fig_path, width = 21, height = 15, units = "cm")

b_new <- bias_sum_plot + theme(axis.text.x = element_blank())
(b_new / mse_sum_plot) + plot_layout(guides = "collect")
ggsave("miss_bin_bias_mse.png", path = fig_path, width = 21, height = 15, units = "cm")

# trying to get all the axes to vary so we can see the variation better?

make_bm_plots <- function(scen) {
  b_plot <- metrics_summary %>%
    filter(scenario == scen) %>%
    ggplot(aes(x = model, y = mean_bias, colour = method, ymin = mean_bias - mcse_bias, ymax = mean_bias + mcse_bias)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbar(position = position_dodge(width = 0.5), linewidth = 0.3) +
    facet_grid(rows = vars(mechanism), cols = vars(scenario), scales = "free") +
    scale_colour_paletteer_d("rcartocolor::Safe") +
    theme_minimal() +
    labs(y = "Bias",
         x = "Model") +
    theme(axis.text.x = element_blank())
  
  m_plot <- metrics_summary %>%
    filter(scenario == scen) %>%
    ggplot(aes(x = model, y = mean_mse, colour = method, ymin = mean_mse - mcse_mse, ymax = mean_mse + mcse_mse)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbar(position = position_dodge(width = 0.5), linewidth = 0.3) +
    facet_grid(rows = vars(mechanism), cols = vars(scenario), scales = "free_y") +
    scale_colour_paletteer_d("rcartocolor::Safe") +
    theme_minimal() +
    labs(y = "MSE",
         x = "Model") +
    theme(axis.text = element_blank())
  list(b_plot = b_plot, m_plot = m_plot)
}

scens <- metrics %>% select(scenario) %>% unique() %>% pull(scenario)

sum_plots <- pmap(list(scens), make_bm_plots)

(sum_plots[[1]]$b_plot + sum_plots[[2]]$b_plot + sum_plots[[3]]$b_plot + sum_plots[[4]]$b_plot) + plot_layout(guides = "collect")
