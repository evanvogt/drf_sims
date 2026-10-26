##########
# title: figures for the thesis - the competing-risk SNR check, all scenarios
##########
# snr_results.qmd's first figure (RMSE ratio, DR / causal forest on PV) for
# every scenario at once. Each point is a cell's geometric mean of the per-run
# ratio, exp(mean log ratio), with a 95% CI from the MCSE of the mean log ratio
# (snr_ratios.RDS, competing_risk/snr/snr_metrics.R); below 1, DR did better.
#
# Two panels per event, as the qmd:
# - "More signal": the n = 500 cells against their population SD(tau) /
#   SD(noise) (snr_population.RDS) on a log scale, one line per scenario
#   through its HTE multipliers, so the scenarios share one x axis. Scenarios
#   1-2 have no HTE multiplier, so they are a single point.
# - "More data": the HTE x1 cells against n, on a log scale.
# The baseline (n = 500, HTE x1) is in both. Colour is the scenario, rows are
# the arms. Scenarios are dodged slightly on x so they don't overlap.
#
# SNR-check results only (100 runs per cell, censoring on, rho = 0.5); the
# main study's figure is rel_eff_snr.R.
#
# Writes to ../results/thesis_figures/snr/:
#   snr_checks.png       DR-RF on PV, DR-RSF and the oracle DR-learner
#   snr_checks_semi.png  their semi-oracle twins (true propensity 0.5)
#
#   cd results_processing/thesis_figures; Rscript snr_checks.R
#   Rscript snr_checks.R [SNR results dir] [figure dir]   # mock smoke run

library(here)
source(here("R", "figures.R"))

# paths
args <- commandArgs(trailingOnly = TRUE)
res_path <- if (length(args) >= 1) {
  args[1]
} else {
  file.path(dirname(here()), "results", "competing_risk_snr")
}
fig_path <- if (length(args) >= 2) {
  args[2]
} else {
  file.path(dirname(here()), "results", "thesis_figures", "snr")
}
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

# the figures' rows, in order
ARM_SETS <- list(
  main = c("pseudo_dr_whole_oob", "rsf_dr_oob", "pseudo_dr_oracle"),
  semi = c("pseudo_dr_whole_oob_semi", "rsf_dr_oob_semi")
)
# snr_results.qmd's labels
ARM_LABELS <- c(
  SURV_ARM_LABELS[c("pseudo_dr_whole_oob", "rsf_dr_oob")],
  pseudo_dr_oracle = "DR-oracle (true μ)",
  pseudo_dr_whole_oob_semi = "DR-RF on PV (e = 0.5)",
  rsf_dr_oob_semi = "DR-RSF (e = 0.5)"
)
PANEL_LABELS <- c(
  signal = "More signal: SD(τ) / SD(noise)\nat n = 500, log scale",
  data = "More data: n at HTE ×1,\nlog scale"
)
Y_LAB <- "RMSE ratio, DR / causal forest on PV\n(geometric mean, 95% CI)"
Z <- qnorm(0.975)
# the scenario dodge, as a share of the panel's x range
X_DODGE <- 0.008

# ---- data --------------------------------------------------------------------

ratios <- readRDS(file.path(res_path, "snr_ratios.RDS"))
population <- readRDS(file.path(res_path, "snr_population.RDS"))
# written before the check had more than scenario 7 (snr_population.R before
# 2026-10-25)
if (!"scenario" %in% names(population)) {
  stop("snr_population.RDS has no scenario column: rerun ",
       "competing_risk/snr/snr_population.R", call. = FALSE)
}

N_LEVELS <- sort(unique(ratios$n))

est <- ratios %>%
  mutate(est = geo_mean,
         lo = exp(log(geo_mean) - Z * mcse_log),
         hi = exp(log(geo_mean) + Z * mcse_log)) %>%
  left_join(select(population, scenario, hte_mult, target, snr),
            by = c("scenario", "hte_mult", "target"))

# each cell on the panel(s) it belongs to, x on the log10 scale; the baseline
# is in both
est <- bind_rows(
  est %>% filter(n == 500) %>% mutate(panel = "signal", x = log10(snr)),
  est %>% filter(hte_mult == 1) %>% mutate(panel = "data", x = log10(n))
) %>%
  group_by(panel) %>%
  mutate(x = x + (scenario - mean(unique(scenario))) * X_DODGE * diff(range(x))) %>%
  ungroup() %>%
  mutate(panel = factor(PANEL_LABELS[panel], levels = PANEL_LABELS),
         scenario = factor(paste0("S", scenario),
                           levels = paste0("S", sort(unique(scenario)))))

missing_snr <- est %>% filter(is.na(x))
if (nrow(missing_snr) > 0) {
  stop(nrow(missing_snr), " cells have no population SNR: rerun ",
       "competing_risk/snr/snr_population.R", call. = FALSE)
}

# ---- figures -----------------------------------------------------------------

# One x scale serves both panels (scales = "free_x"), and each panel's limits
# say which it is: log10(n) is above 2, log10(SD(tau) / SD(noise)) is below it
x_breaks <- function(lims) {
  if (lims[1] > 2) return(log10(N_LEVELS))
  log10(scales::breaks_log(n = 6)(10^lims))
}
x_labels <- function(b) {
  vapply(b, function(v) {
    if (is.na(v)) return(NA_character_)
    if (v > 2) format(round(10^v), big.mark = ",") else
      format(signif(10^v, 2), drop0trailing = TRUE, scientific = FALSE)
  }, character(1))
}

#' The ratio figure for one set of arms: a row per arm, Event 1 and Event 2
#' side by side, each with its two panels
snr_checks_plot <- function(arms) {
  d <- est %>%
    filter(framework %in% arms) %>%
    mutate(arm = factor(ARM_LABELS[framework], levels = ARM_LABELS[arms]))
  # each scenario's label at its right-hand end of the "more signal" panel
  ends <- d %>%
    filter(panel == PANEL_LABELS[["signal"]]) %>%
    group_by(arm, target, panel, scenario) %>%
    slice_max(x, n = 1) %>%
    ungroup()

  ggplot(d, aes(x, est, colour = scenario, group = scenario)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_line(linewidth = 0.4) +
    geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.5) +
    geom_point(size = 1.5) +
    geom_text(data = ends, aes(label = scenario), hjust = -0.4, size = 2.3,
              show.legend = FALSE) +
    facet_grid(arm ~ target + panel, scales = "free_x", space = "free_x") +
    scale_x_continuous(breaks = x_breaks, labels = x_labels,
                       minor_breaks = NULL,
                       expand = expansion(mult = c(0.05, 0.12))) +
    scale_y_log10() +
    drf_scale() +
    labs(x = NULL, y = Y_LAB, colour = "Scenario") +
    drf_theme() +
    theme(legend.position = "bottom",
          panel.grid.major.x = element_blank()) +
    guides(colour = guide_legend(nrow = 1))
}

figs <- list()
for (set in names(ARM_SETS)) {
  if (!any(ARM_SETS[[set]] %in% est$framework)) {
    message("no ", set, " arms in snr_ratios.RDS - figure skipped")
    next
  }
  figs[[set]] <- snr_checks_plot(ARM_SETS[[set]])
  name <- paste0("snr_checks", if (set != "main") paste0("_", set), ".png")
  save_fig(name, fig_path, width = 28,
           height = 5 + 5 * length(ARM_SETS[[set]]), plot = figs[[set]])
}
