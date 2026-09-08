##########
# title: nuisance / pseudo-value diagnostic figures - competing risk
##########
# Reads the two objects surv_nuisance_extract.R produces and plots them:
#   - propensity overlap (W.hat by treatment) - the standard positivity check
#   - the estimated pseudo-value regressions (pseudo.hat/pseudo0.hat/pseudo.hat.cf)
#   - the DR pseudo-outcome (po) distribution, which the AIPW correction term
#     can make heavy-tailed when W.hat sits near the trim bounds
#   - a po-vs-W.hat scatter, to see whether the two line up
#   - a full-corpus (all runs, not just the sampled ones) overview of the
#     propensity trimming rate and the po extreme-value rate
#
# Facets are already crowded (7 scenarios x 5 arms x 3 estimands x 2
# censoring), so each plot below picks one slice of that (one censoring value,
# sometimes one estimand) rather than crossing everything into one panel.

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
source(here("R", "figures.R"))
source(here("competing_risk/surv_config.R"))

# rcartocolor::Safe as a fill scale, to match drf_scale()'s colour version -
# not added to R/figures.R since every other helper there uses colour, not
# fill (distribution_plot() etc. use fill = "transparent" on purpose).
nuisance_fill_scale <- function() scale_fill_paletteer_d("rcartocolor::Safe")

NUISANCE_ARM_LABELS <- c(
  rf_whole_oob = "RF (whole, OOB)",
  rf_whole_scf = "RF (whole, single crossfit)",
  rf_cvps_scf = "RF (cvps, single crossfit)",
  sl_whole = "SuperLearner (whole)",
  sl_cvps = "SuperLearner (cvps)"
)

fig_dir <- file.path(study$res_path, "nuisance_figs")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

indiv <- readRDS(file.path(study$res_path, "nuisance_indiv_sample.RDS")) %>%
  mutate(
    scenario = factor(scenario),
    arm_label = label_factor(arm_label, NUISANCE_ARM_LABELS),
    estimand = factor(estimand, levels = c("RMTL1", "RMTL2", "RMSTc"))
  )

run_summary <- readRDS(file.path(study$res_path, "nuisance_run_summary.RDS")) %>%
  mutate(
    scenario = factor(scenario),
    arm_label = label_factor(arm_label, NUISANCE_ARM_LABELS),
    estimand = factor(estimand, levels = c("RMTL1", "RMTL2", "RMSTc"))
  )

# ---- plot builders -----------------------------------------------------------

#' Propensity overlap: density of W.hat by actual treatment
#'
#' W.hat is refit per estimand (surv_models.R's by_estimand()) but its shape
#' barely changes across the three - RMTL1 is used here as the representative
#' one so scenario x arm stays a single, readable grid.
propensity_overlap_plot <- function(indiv, censor_value) {
  plot_df <- filter(indiv, censoring == censor_value, estimand == "RMTL1")

  ggplot(plot_df, aes(x = W.hat, fill = factor(W))) +
    geom_density(alpha = 0.5, colour = NA) +
    geom_vline(xintercept = c(0.05, 0.95), linetype = "dashed", colour = "grey40") +
    facet_grid(scenario ~ arm_label) +
    nuisance_fill_scale() +
    drf_theme() +
    labs(
      title = paste0("Propensity score overlap (censoring = ", censor_value, ")"),
      x = "W.hat", y = "Density", fill = "Treatment (W)"
    )
}

#' Estimated pseudo-value regressions, by nuisance arm
pseudo_value_plot <- function(indiv, censor_value, estimand_value) {
  plot_df <- indiv %>%
    filter(censoring == censor_value, estimand == estimand_value) %>%
    pivot_longer(
      c(pseudo.hat, pseudo0.hat, pseudo.hat.cf),
      names_to = "quantity",
      values_to = "value"
    )

  ggplot(plot_df, aes(x = arm_label, y = value, colour = base_learner)) +
    geom_boxplot(fill = "transparent", outlier.shape = NA) +
    facet_grid(scenario ~ quantity, scales = "free_y") +
    drf_scale() +
    drf_theme(rotate_x = TRUE) +
    labs(
      title = paste0("Estimated pseudo-values (", estimand_value, ", censoring = ", censor_value, ")"),
      x = "Nuisance arm", y = estimand_value, colour = "Base learner"
    )
}

#' DR pseudo-outcome (po) distribution, by nuisance arm
po_distribution_plot <- function(indiv, censor_value) {
  plot_df <- filter(indiv, censoring == censor_value)

  ggplot(plot_df, aes(x = arm_label, y = po, colour = arm_label)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_boxplot(fill = "transparent", outlier.shape = NA) +
    facet_grid(scenario ~ estimand, scales = "free_y") +
    drf_scale() +
    drf_theme(rotate_x = TRUE) +
    labs(
      title = paste0("DR pseudo-outcome (po) distribution (censoring = ", censor_value, ")"),
      x = "Nuisance arm", y = "po"
    ) +
    theme(legend.position = "none")
}

#' po vs W.hat - do extreme pseudo-outcomes coincide with extreme propensities?
po_vs_what_plot <- function(indiv, censor_value, estimand_value) {
  plot_df <- filter(indiv, censoring == censor_value, estimand == estimand_value)

  ggplot(plot_df, aes(x = W.hat, y = po)) +
    geom_bin2d(bins = 40) +
    geom_vline(xintercept = c(0.05, 0.95), linetype = "dashed", colour = "white") +
    facet_grid(scenario ~ arm_label) +
    scale_fill_viridis_c() +
    drf_theme() +
    labs(
      title = paste0("po vs W.hat (", estimand_value, ", censoring = ", censor_value, ")"),
      x = "W.hat", y = "po", fill = "Count"
    )
}

#' Propensity trimming rate across every run (not just the sampled ones)
w_hat_trim_plot <- function(run_summary) {
  plot_df <- run_summary %>%
    filter(quantity == "W.hat") %>%
    mutate(trim_rate = prop_at_lo + prop_at_hi)

  ggplot(plot_df, aes(x = arm_label, y = trim_rate, colour = arm_label)) +
    geom_boxplot(fill = "transparent", outlier.shape = NA) +
    facet_grid(censoring ~ scenario) +
    scale_y_continuous(labels = percent) +
    drf_scale() +
    drf_theme(rotate_x = TRUE) +
    labs(
      title = "Propensity trimming rate across all runs",
      x = "Nuisance arm", y = "Share of individuals at the 0.05/0.95 trim bound"
    ) +
    theme(legend.position = "none")
}

#' Share of extreme po values across every run (not just the sampled ones)
po_extreme_rate_plot <- function(run_summary) {
  plot_df <- filter(run_summary, quantity == "po")

  ggplot(plot_df, aes(x = arm_label, y = prop_extreme, colour = arm_label)) +
    geom_boxplot(fill = "transparent", outlier.shape = NA) +
    facet_grid(censoring ~ scenario) +
    scale_y_continuous(labels = percent) +
    drf_scale() +
    drf_theme(rotate_x = TRUE) +
    labs(
      title = "Share of extreme DR pseudo-outcome (po) values across all runs",
      x = "Nuisance arm", y = "Share with |po - median| > 3xIQR"
    ) +
    theme(legend.position = "none")
}

# ---- render and save ----------------------------------------------------------

for (cv in c(TRUE, FALSE)) {
  cv_label <- if (cv) "censored" else "uncensored"

  save_fig(
    paste0("propensity_overlap_", cv_label, ".png"),
    path = fig_dir, width = 30, height = 20,
    plot = propensity_overlap_plot(indiv, cv)
  )

  save_fig(
    paste0("po_distribution_", cv_label, ".png"),
    path = fig_dir, width = 30, height = 20,
    plot = po_distribution_plot(indiv, cv)
  )

  for (est in c("RMTL1", "RMTL2", "RMSTc")) {
    save_fig(
      paste0("pseudo_values_", est, "_", cv_label, ".png"),
      path = fig_dir, width = 30, height = 20,
      plot = pseudo_value_plot(indiv, cv, est)
    )
    save_fig(
      paste0("po_vs_what_", est, "_", cv_label, ".png"),
      path = fig_dir, width = 30, height = 20,
      plot = po_vs_what_plot(indiv, cv, est)
    )
  }
}

save_fig(
  "w_hat_trim_rate_overview.png",
  path = fig_dir, width = 30, height = 15,
  plot = w_hat_trim_plot(run_summary)
)
save_fig(
  "po_extreme_rate_overview.png",
  path = fig_dir, width = 30, height = 15,
  plot = po_extreme_rate_plot(run_summary)
)

message("Nuisance figures written to ", fig_dir)
