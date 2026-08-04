##########
# title: figures for the thesis chapter - competing risk
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
res_path <- file.path(dirname(path), "results", "competing_risk")
fig_path <- file.path(dirname(path), "results", "thesis_figures", "surv")
dir.create(fig_path, showWarnings = F, recursive = T)

# data
metrics <- readRDS(file.path(res_path, "surv_metrics.RDS"))

# clean up
metrics <- metrics %>%
  mutate(
    description = factor(
      case_match(scenario,
                 1 ~ "ATE on EOI only",
                 2 ~ "ATE on CE only",
                 3 ~ "HTE on EOI only",
                 4 ~ "HTE on EOI, ATE on CE",
                 5 ~ "HTE on CE only",
                 6 ~ "HTE on CE, ATE on EOI",
                 7 ~ "HTE on both"),
      levels = c("ATE on EOI only", "ATE on CE only", "HTE on EOI only",
                 "HTE on EOI, ATE on CE", "HTE on CE only",
                 "HTE on CE, ATE on EOI", "HTE on both")),
    scenario = factor(scenario, levels = c(1:7)),
    framework = factor(
      case_match(framework,
                 "csf_cs" ~ "CSF - censoring CEs",
                 "csf_sh" ~ "CSF - subdistribution",
                 "ipw"      ~ "IPW",
                 "pseudo_cf" ~ "pseudovalue CF",
                 "pseudo_dr" ~ "pseudovalue DR RF"),
      levels = c("CSF - censoring CEs", "CSF - subdistribution", "IPW", "pseudovalue CF", "pseudovalue DR RF")),
    target = factor(target, levels = c("Event 1", "Event 2", "Combined"))
  )

# per scenario summaries
metrics_summary <- metrics %>%
  group_by(scenario, description, n, censoring, framework, target) %>%
  summarise(
    mean_bias = mean(bias, na.rm = T),
    mcse_bias = sd(bias, na.rm = T)/sqrt(n()),
    mean_mse = mean(mse, na.rm = T),
    mcse_mse = sd(mse, na.rm = T)/sqrt(n()),
    mean_corr = mean(corr, na.rm = T),
    mcse_corr = sd(corr, na.rm = T)/sqrt(n()),
    .groups = "drop"
  )

# create a cleaner metrics table to save and look at to write results
metrics_sum_tidy <- metrics_summary %>%
  mutate(
    bias_lb = signif(mean_bias - mcse_bias, 4),
    bias_ub = signif(mean_bias + mcse_bias, 4),
    mse_lb = signif(mean_mse - mcse_mse, 4),
    mse_ub = signif(mean_mse + mcse_mse, 4),
    corr_lb = signif(mean_corr - mcse_corr, 4),
    corr_ub = signif(mean_corr + mcse_corr, 4),
    bias = signif(mean_bias, 4),
    mse = signif(mean_mse, 4),
    corr = signif(mean_corr, 4)
  ) %>%
  mutate(
    bias_tidy = paste0(bias, " (", bias_lb, ", ", bias_ub, ")"),
    mse_tidy = paste0(mse, " (", mse_lb, ", ", mse_ub, ")"),
    corr_tidy = paste0(corr, " (", corr_lb, ", ", corr_ub, ")")
  ) %>%
  select(scenario, description, censoring, target, framework, bias, bias_lb, bias_ub, bias_tidy, mse, mse_lb, mse_ub, mse_tidy, corr, corr_lb, corr_ub, corr_tidy)
write.csv(metrics_sum_tidy, file.path(res_path, "surv_metric_sum_tidy.csv"))

bias_plot <- metrics_summary %>%
  filter(!(framework == "IPW" & target == "Event 2")) %>%
  #filter(!(framework %in% c("CSF - censoring CEs", "IPW"))) %>%
  ggplot(aes(x = censoring, colour = framework, y = mean_bias, ymin = mean_bias - mcse_bias, ymax = mean_bias + mcse_bias)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(position = position_dodge(width = 0.5), linewidth = 0.3) +
  facet_grid(cols = vars(scenario), rows = vars(target), scales = "free_y") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "Bias",
       x = "Censoring") +
  theme(axis.text.x = element_text(angle = 90))
ggsave("surv_bias.png", path = fig_path, width = 21, height = 15, units = "cm")

mse_plot <- metrics_summary %>%
  filter(!(framework == "IPW" & target == "Event 2")) %>%
  #filter(!(framework %in% c("CSF - censoring CEs", "IPW"))) %>%
  ggplot(aes(x = censoring, colour = framework, y = mean_mse, ymin = mean_mse - mcse_mse, ymax = mean_mse + mcse_mse)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(position = position_dodge(width = 0.5), linewidth = 0.3) +
  facet_grid(cols = vars(scenario), rows = vars(target), scales = "free_y") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "MSE",
       x = "Censoring") +
  theme(axis.text.x = element_text(angle = 90))
ggsave("surv_mse.png", path = fig_path, width = 21, height = 15, units = "cm")


# scratch
metrics_summary %>%
  filter(framework == "IPW") %>%
  #filter(framework %in% c("CSF - censoring CEs", "CSF - subdistribution")) %>%
  ggplot(aes(x = censoring, colour = framework, y = mean_bias, ymin = mean_bias - mcse_bias, ymax = mean_bias + mcse_bias)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(position = position_dodge(width = 0.5), linewidth = 0.3) +
  facet_grid(cols = vars(scenario), rows = vars(target), scales = "free_y") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "Bias",
       x = "Censoring") +
  theme(axis.text.x = element_text(angle = 90))

metrics_summary %>%
  filter(framework == "IPW") %>%
  #filter(framework %in% c("CSF - censoring CEs", "CSF - subdistribution")) %>%
  ggplot(aes(x = censoring, colour = framework, y = mean_mse, ymin = mean_mse - mcse_mse, ymax = mean_mse + mcse_mse)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(position = position_dodge(width = 0.5), linewidth = 0.3) +
  facet_grid(cols = vars(scenario), rows = vars(target), scales = "free_y") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "MSE",
       x = "Censoring") +
  theme(axis.text.x = element_text(angle = 90))
