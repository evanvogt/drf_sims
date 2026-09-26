###############
# ICTMC presentation comprisk figs
##############

library(here)
library(tidyverse)
library(patchwork)
library(ggridges)
source(here("R", "figures.R"))

path <- here()
figdir <- file.path(dirname(path), "results", "ICTMC_figs", "presentation")
dir.create(figdir, showWarnings = F)


metrics <- readRDS(here("../collected_metrics/competing_risk/surv_metrics.RDS"))

# labels -----
scenario_labels <- c(
  `1` = "Constant CATE Event 1",
  `3` = "Varying CATE Event 1",
  `4` = "Varying CATE Event 1,\nConstant CATE Event 2",
  `2` = "Constant CATE Event 2",
  `5` = "Varying CATE Event 2",
  `6` = "Varying CATE Event 2,\nConstant CATE Event 1",
  `7` = "Varying CATEs both events"
)
framework_family <- c(
  ipw = "IPW",
  csf_cs = "CSF (cause-specific)",
  csf_sh = "Adapted CSF",
  pseudo_cf_whole_oob = "Pseudovalue CF",
  pseudo_cf_whole_scf = "Pseudovalue CF",
  pseudo_cf_cvps_scf = "Pseudovalue CF",
  pseudo_dr_whole_oob = "Pseudovalue DR RF",
  pseudo_dr_whole_scf = "Pseudovalue DR RF",
  pseudo_dr_cvps_scf = "Pseudovalue DR RF",
  sl_t_whole = "SL T-learner",
  sl_t_cvps = "SL T-learner",
  sl_t_split = "SL T-learner",
  sl_dr_whole = "Pseudovalue DR SL",
  sl_dr_cvps = "Pseudovalue DR SL"
)
family_levels <- c(
  "IPW",
  "CSF (cause-specific)",
  "Adapted CSF",
  "Pseudovalue CF",
  "Pseudovalue DR RF",
  "SL T-learner",
  "Pseudovalue DR SL"
)

framework_variant <- c(
  ipw = "single",
  csf_cs = "single",
  csf_sh = "single",
  pseudo_cf_whole_oob = "whole_oob",
  pseudo_cf_whole_scf = "whole_scf",
  pseudo_cf_cvps_scf = "cvps_scf",
  pseudo_dr_whole_oob = "whole_oob",
  pseudo_dr_whole_scf = "whole_scf",
  pseudo_dr_cvps_scf = "cvps_scf",
  sl_t_whole = "whole",
  sl_t_cvps = "cvps",
  sl_t_split = "split",
  sl_dr_whole = "whole",
  sl_dr_cvps = "cvps"
)
variant_levels <- c(
  "single",
  "whole_oob",
  "whole_scf",
  "cvps_scf",
  "whole",
  "cvps",
  "split"
)
variant_labels <- c(
  single = "(single arm)",
  whole_oob = "Whole-sample, OOB",
  whole_scf = "Whole-sample, single CF",
  cvps_scf = "Crossfit PV, single CF",
  whole = "Whole-sample PV",
  cvps = "Crossfit PV",
  split = "Split PV"
)

# colour palette - same pair as ictmc_figs.R's threshold_palette (subgroup
# validation significance thresholds), reused here for the censoring/
# no-censoring split so the two talks share a colour vocabulary
censoring_palette <- c(
  "Censoring" = "#dc143c",
  "No Censoring" = "#4b0082"
)

# tidy metrics -----
metrics <- metrics %>%
  mutate(
    family = factor(
      recode(framework, !!!framework_family),
      levels = family_levels
    ),
    variant = recode(framework, !!!framework_variant),
    variant = factor(
      recode(variant, !!!variant_labels),
      levels = unname(variant_labels[variant_levels])
    ),
    scenario = factor(
      scenario,
      levels = names(scenario_labels) %>% as.integer(),
      labels = unname(scenario_labels)
    ),
    censoring = factor(
      censoring,
      levels = c(TRUE, FALSE),
      labels = c("Censoring", "No Censoring")
    ),
    target = factor(target, levels = c("Event 1", "Event 2", "Combined"))
  ) %>%
  droplevels()

# filter metrics -----
metrics <- metrics %>%
  filter(
    !(family %in%
      c("IPW", "CSF (cause-specific)", "SL T-learner", "Pseudovalue DR SL"))
  ) %>%
  filter(target != "Combined") %>%
  filter(
    variant %in% c("(single arm)", "Crossfit PV, single CF", "Crossfit PV")
  ) %>%
  droplevels()

metrics_summary <- metrics %>%
  group_by(scenario, censoring, target, family, variant) %>%
  summarise(
    mean_bias = mean(bias, na.rm = T),
    mcse_bias = sd(bias, na.rm = T) / sqrt(sum(!is.na(bias))),
    mean_ate_bias = mean(ate_bias, na.rm = T),
    mcse_ate_bias = sd(ate_bias, na.rm = T) / sqrt(sum(!is.na(ate_bias))),
    mean_rel_ate_bias = mean(rel_ate_bias, na.rm = T),
    mcse_rel_ate_bias = sd(rel_ate_bias, na.rm = T) /
      sqrt(sum(!is.na(rel_ate_bias))),
    mean_rel_bias_cate = mean(rel_bias_cate, na.rm = T),
    mcse_rel_bias_cate = sd(rel_bias_cate, na.rm = T) /
      sqrt(sum(!is.na(rel_bias_cate))),
    mean_mse = mean(mse, na.rm = T),
    mcse_mse = sd(mse, na.rm = T) / sqrt(sum(!is.na(mse))),
    mean_rmse = mean(rmse, na.rm = T),
    mcse_rmse = sd(rmse, na.rm = T) / sqrt(sum(!is.na(rmse))),
    mean_mae = mean(mae, na.rm = T),
    mcse_mae = sd(mae, na.rm = T) / sqrt(sum(!is.na(mae))),
    mean_corr = mean(corr, na.rm = T),
    mcse_corr = sd(corr, na.rm = T) / sqrt(sum(!is.na(corr))),
    mean_spearman = mean(spearman, na.rm = T),
    mcse_spearman = sd(spearman, na.rm = T) / sqrt(sum(!is.na(spearman))),
    mean_sign_acc = mean(sign_acc, na.rm = T),
    mcse_sign_acc = sd(sign_acc, na.rm = T) / sqrt(sum(!is.na(sign_acc))),
    mean_c_stat = mean(c_stat, na.rm = T),
    mcse_c_stat = sd(c_stat, na.rm = T) / sqrt(sum(!is.na(c_stat))),
    mean_n_na = mean(n_na, na.rm = T),
    total_n_na = sum(n_na, na.rm = T),
    .groups = "drop"
  )

summary_plot <- function(
  df,
  mean_col,
  mcse_col,
  title,
  ylab,
  hline = 0,
  alpha = 0.05
) {
  z <- qnorm(1 - alpha / 2)
  df %>%
    mutate(
      est = .data[[mean_col]],
      lo = .data[[mean_col]] - z * .data[[mcse_col]],
      hi = .data[[mean_col]] + z * .data[[mcse_col]]
    ) %>%
    ggplot(aes(
      x = family,
      y = est,
      ymin = lo,
      ymax = hi,
      colour = censoring
    )) +
    geom_hline(yintercept = hline, linetype = "dashed") +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbar(
      position = position_dodge(width = 0.5),
      linewidth = 0.3,
      width = 0.3
    ) +
    facet_grid(rows = vars(target), cols = vars(scenario), scales = "free") +
    scale_colour_manual(values = censoring_palette) +
    theme_light() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(colour = "black")
    ) +
    labs(
      title = title,
      y = ylab,
      x = "Model",
      colour = "Censoring"
    )
}

summary_plot(metrics_summary, "mean_bias", "mcse_bias", "", "Mean bias")
ggsave(
  "bias_both_events.png",
  path = figdir,
  height = 20,
  width = 30,
  units = "cm"
)
summary_plot(metrics_summary, "mean_mse", "mcse_mse", "", "Mean MSE")
ggsave(
  "mse_both_events.png",
  path = figdir,
  height = 20,
  width = 30,
  units = "cm"
)

# Event 1 only: bias top row / MSE bottom row -----------
metrics_summary_e1 <- metrics_summary %>%
  filter(target == "Event 1") %>%
  droplevels()

bias_e1_full <- summary_plot(
  metrics_summary_e1,
  "mean_bias",
  "mcse_bias",
  "",
  "Mean bias"
)
ggsave(
  "bias_event1.png",
  bias_e1_full,
  path = figdir,
  height = 10,
  width = 30,
  units = "cm"
)

mse_e1_panel <- summary_plot(
  metrics_summary_e1,
  "mean_mse",
  "mcse_mse",
  "",
  "Mean MSE"
)
ggsave(
  "mse_event1.png",
  mse_e1_panel,
  path = figdir,
  height = 10,
  width = 30,
  units = "cm"
)

# blanked-x variant of the bias panel, just for the stacked bias/MSE combo
# below - the standalone bias_event1.png above keeps its own x-axis
bias_e1_panel <- bias_e1_full +
  labs(x = NULL) +
  theme(axis.text.x = element_blank())

(bias_e1_panel / mse_e1_panel) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(
  "bias_mse_event1.png",
  path = figdir,
  height = 20,
  width = 30,
  units = "cm"
)

# Event 2 only: bias top row / MSE bottom row -----------
metrics_summary_e2 <- metrics_summary %>%
  filter(target == "Event 2") %>%
  droplevels()

bias_e2_full <- summary_plot(
  metrics_summary_e2,
  "mean_bias",
  "mcse_bias",
  "",
  "Mean bias"
)
ggsave(
  "bias_event2.png",
  bias_e2_full,
  path = figdir,
  height = 10,
  width = 30,
  units = "cm"
)

mse_e2_panel <- summary_plot(
  metrics_summary_e2,
  "mean_mse",
  "mcse_mse",
  "",
  "Mean MSE"
)
ggsave(
  "mse_event2.png",
  mse_e2_panel,
  path = figdir,
  height = 10,
  width = 30,
  units = "cm"
)

# blanked-x variant of the bias panel, just for the stacked bias/MSE combo
# below - the standalone bias_event2.png above keeps its own x-axis
bias_e2_panel <- bias_e2_full +
  labs(x = NULL) +
  theme(axis.text.x = element_blank())

(bias_e2_panel / mse_e2_panel) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(
  "bias_mse_event2.png",
  path = figdir,
  height = 20,
  width = 30,
  units = "cm"
)

# Single-event plots split by scenario group -----------
# scenarios 1, 3, 4 vary/hold the Event 1 CATE; 2, 5, 6, 7 vary/hold the
# Event 2 CATE (plus scenario 7, where both vary). All 7 scenario columns in
# one panel is too wide to read on a slide, so split each single-event
# bias/MSE plot into these two narrower panels instead.
scenario_group_e1 <- unname(scenario_labels[c("1", "3", "4")])
scenario_group_e2 <- unname(scenario_labels[c("2", "5", "6", "7")])

single_event_scenario_plots <- function(
  metrics_summary_df,
  target_name,
  scenario_group,
  file_suffix,
  width
) {
  df <- metrics_summary_df %>%
    filter(target == target_name, scenario %in% scenario_group) %>%
    droplevels()
  target_slug <- tolower(gsub(" ", "", target_name))

  bias_panel <- summary_plot(df, "mean_bias", "mcse_bias", "", "Mean bias")
  ggsave(
    paste0("bias_", target_slug, "_", file_suffix, ".png"),
    bias_panel,
    path = figdir,
    height = 10,
    width = width,
    units = "cm"
  )

  mse_panel <- summary_plot(df, "mean_mse", "mcse_mse", "", "Mean MSE")
  ggsave(
    paste0("mse_", target_slug, "_", file_suffix, ".png"),
    mse_panel,
    path = figdir,
    height = 10,
    width = width,
    units = "cm"
  )

  bias_panel_blanked <- bias_panel +
    labs(x = NULL) +
    theme(axis.text.x = element_blank())
  stacked_panel <- (bias_panel_blanked / mse_panel) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  ggsave(
    paste0("bias_mse_", target_slug, "_", file_suffix, ".png"),
    stacked_panel,
    path = figdir,
    height = 20,
    width = width,
    units = "cm"
  )
}

single_event_scenario_plots(
  metrics_summary,
  "Event 1",
  scenario_group_e1,
  "scenarios_1_3_4",
  width = 17
)
single_event_scenario_plots(
  metrics_summary,
  "Event 1",
  scenario_group_e2,
  "scenarios_2_5_6_7",
  width = 23
)
single_event_scenario_plots(
  metrics_summary,
  "Event 2",
  scenario_group_e1,
  "scenarios_1_3_4",
  width = 17
)
single_event_scenario_plots(
  metrics_summary,
  "Event 2",
  scenario_group_e2,
  "scenarios_2_5_6_7",
  width = 23
)
