##########
# title: figures for the thesis chapter - sample size, correlated covariates
##########
# The main chapter's sample-size figures: the correlated-covariate studies
# (sample_size/correlated/) at rho = 0.5, scenarios 1-4. The same figures for
# scenarios 5-10, also at rho = 0.5, are the appendix's (*_corr_supp_*.png;
# the appendix's tables, at both rhos, are thesis_tables/ss_tables.R and
# ss_test_tables.R). The independent-covariate studies (cts_ss.R, bin_ss.R)
# are the appendix's too.
#
# Two figures per outcome and scenario set, the scenarios across in each:
#   - CATE bias on the top row, RMSE on the bottom;
#   - each HTE test's rejection rate at 0.1, one row per test, with the same
#     tests run on the true values as a reference series.
# Labels, palette, summaries and figure sizing come from R/figures.R. This
# script carries only the paths and this study's filters.

library(here)
library(patchwork)
source(here("R", "figures.R"))

# paths
path <- here()
res_path <- file.path(dirname(path), "results", "correlated")
fig_path <- file.path(dirname(path), "results", "thesis_figures", "sample_size")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

RHO <- 0.5

# the appendix's scenarios: six panels across rather than four, so wider
SUPP_SCENARIOS <- 5:10
SUPP_WIDTH <- 29

#' Bias (top) and RMSE (bottom) by sample size, one panel per scenario
#'
#' @param outcome results folder, "continuous" or "binary"
#' @param prefix file prefix, "cts" or "bin"
#' @param scale_lab appended to the y labels, e.g. " (risk difference)"
#' @param scenarios,labels the scenarios to plot, and their strip labels
#' @param suffix inserted into the file name, e.g. "_supp"
#' @param width figure width in cm, as save_fig()
bias_rmse_figure <- function(outcome, prefix, scale_lab = "",
                             scenarios = 1:4, labels = SS_SCENARIO_LABELS,
                             suffix = "", width = 21) {
  metrics <- readRDS(file.path(res_path, outcome,
                               paste0(prefix, "_corr_metrics.RDS"))) %>%
    filter(rho == RHO, scenario %in% scenarios) %>%
    apply_labels(labels) %>%
    mutate(n = factor(n, levels = c(100, 250, 500, 1000)))

  # rmse isn't in summarise_metrics()'s default cols
  metrics_summary <- summarise_metrics(metrics, c("scenario", "n", "model"),
                                       cols = c(bias = "bias", rmse = "rmse"))

  bias_row <- point_range_plot(metrics_summary, "bias",
                               paste0("Bias of the CATE", scale_lab),
                               x = "n", colour = "model", facet_rows = NULL,
                               facet_cols = "scenario", line = TRUE, ci_alpha = 0.7, hline = 0,
                               blank_x = TRUE) +
    labs(x = NULL)

  # scenario names are on the top row's strips already
  rmse_row <- point_range_plot(metrics_summary, "rmse",
                               paste0("RMSE of the CATE", scale_lab),
                               x = "n", colour = "model", facet_rows = NULL,
                               facet_cols = "scenario", line = TRUE, ci_alpha = 0.7,
                               hline = NULL) +
    labs(x = "Sample size") +
    theme(strip.text.x = element_blank(),
          strip.background.x = element_blank())

  fig <- (bias_row / rmse_row) +
    plot_layout(guides = "collect") &
    labs(colour = "Model") &
    theme(legend.position = "bottom") &
    guides(colour = guide_legend(nrow = 2))
  save_fig(paste0(prefix, "_corr", suffix, "_bias_rmse.png"), fig_path,
           width = width, plot = fig)
  fig
}

cts_fig <- bias_rmse_figure("continuous", "cts")
bin_fig <- bias_rmse_figure("binary", "bin")

# the HTE tests' significance threshold, and the power each should reach where
# there is HTE (a conventional target - nothing in the DGM aims for it)
ALPHA <- 0.1
POWER_TARGET <- 0.9

# the HTE tests, in facet-row order
HTE_TESTS <- c(
  BLP_p_os = "BLP (one-sided, HC3)",
  BLP_p = "BLP (two-sided)",
  indep_po = "PO independence",
  indep_cate = "CATE independence"
)

TRUE_LAB <- "True values"

#' Rejection rate at ALPHA by sample size, one row per test, one column per
#' scenario
#'
#' Dashed reference line at ALPHA in the null scenario (nominal size) and at
#' POWER_TARGET in the others.
#'
#' The "True values" series is the BLP and CATE independence tests run on the
#' true CATE and true nuisances (<prefix>_corr_true_cate_tests.RDS, see
#' true_cate_test_row() in R/cate_models.R). The true CATE is constant in the
#' null scenario, so those tests have no true series there. There is no
#' true-value PO test: DR-oracle's own indep_po already tests the true
#' pseudo-outcome, so it stands in as the PO row's true series rather than as a
#' model of its own. Causal forest and the T-learners have no indep_po (see
#' hte_test_metrics() in R/metrics.R), so they draw nothing in that row.
#'
#' @param outcome results folder, "continuous" or "binary"
#' @param prefix file prefix, "cts" or "bin"
#' @param scenarios,labels,suffix,width as bias_rmse_figure()
rejection_figure <- function(outcome, prefix,
                             scenarios = 1:4, labels = SS_SCENARIO_LABELS,
                             suffix = "", width = 21) {
  read_corr <- function(stem) {
    readRDS(file.path(res_path, outcome, paste0(prefix, "_corr_", stem, ".RDS"))) %>%
      filter(rho == RHO, scenario %in% scenarios)
  }
  metrics <- read_corr("metrics")
  true_tests <- read_corr("true_cate_tests")

  est <- metrics %>%
    transmute(scenario, n, model = as.character(model),
              across(all_of(names(HTE_TESTS)))) %>%
    pivot_longer(all_of(names(HTE_TESTS)), names_to = "test", values_to = "p") %>%
    filter(!(model == "dr_oracle" & test == "indep_po"))

  true_cols <- c("BLP_p", "BLP_p_os", "indep_cate")
  truth <- bind_rows(
    true_tests %>%
      select(scenario, n, all_of(true_cols)) %>%
      pivot_longer(all_of(true_cols), names_to = "test", values_to = "p"),
    metrics %>%
      filter(model == "dr_oracle") %>%
      transmute(scenario, n, test = "indep_po", p = indep_po)
  ) %>%
    mutate(model = TRUE_LAB)

  # NA p-values (a constant CATE, or a test with no value for that model) are
  # left out, as summarise_metrics()'s na.rm would anyway
  rejections <- bind_rows(est, truth) %>%
    filter(!is.na(p)) %>%
    mutate(rej = as.numeric(p < ALPHA),
           test = factor(test, levels = names(HTE_TESTS), labels = HTE_TESTS)) %>%
    apply_labels(labels) %>%
    mutate(n = factor(n, levels = c(100, 250, 500, 1000)))

  # rejection is a per-run 0/1 indicator, so the binomial MCSE
  rej_summary <- summarise_metrics(rejections, c("test", "scenario", "n", "model"),
                                   cols = c(rej = "rej"), binomial = "rej",
                                   count_na = character())

  # the models keep their colours from the bias/RMSE figure (Safe, in level
  # order); the true series, label_factor()'s last level, is black
  model_levels <- levels(rej_summary$model)
  pal <- scale_colour_manual(values = setNames(
    c(as.character(paletteer_d("rcartocolor::Safe", length(model_levels) - 1)),
      "black"),
    model_levels
  ))

  # nominal size where there is no HTE, the power target everywhere else
  # (scenarios 5-10 all have HTE)
  scenario_levels <- levels(rej_summary$scenario)
  ref_lines <- tibble(
    scenario = factor(scenario_levels, levels = scenario_levels),
    yintercept = if_else(scenario_levels == SS_SCENARIO_LABELS[["1"]],
                         ALPHA, POWER_TARGET)
  )

  fig <- point_range_plot(rej_summary, "rej", paste("Rejection rate at", ALPHA),
                          x = "n", colour = "model", facet_rows = "test",
                          facet_cols = "scenario", facet_scales = "fixed",
                          line = TRUE, ci_alpha = 0.7, hline = ref_lines,
                          palette = pal) +
    labs(x = "Sample size", colour = "Model") +
    theme(legend.position = "right") #+ # figures are too long with legend at the bottom
    #guides(colour = guide_legend(nrow = 2))
  save_fig(paste0(prefix, "_corr", suffix, "_rejection.png"), fig_path,
           width = width, height = 18, plot = fig)
  fig
}

cts_rej_fig <- rejection_figure("continuous", "cts")
bin_rej_fig <- rejection_figure("binary", "bin")

# ---- appendix: scenarios 5-10 --------------------------------------------------
# Last, and per outcome only once its metrics have scenarios 5-10 (they were
# added to the study after 1-4 had run), so the main figures above never wait
# on them.

for (o in list(c(outcome = "continuous", prefix = "cts"),
               c(outcome = "binary", prefix = "bin"))) {
  has_supp <- readRDS(file.path(res_path, o[["outcome"]],
                                paste0(o[["prefix"]], "_corr_metrics.RDS"))) %>%
    filter(rho == RHO, scenario %in% SUPP_SCENARIOS) %>%
    nrow() > 0
  if (!has_supp) {
    message(o[["outcome"]], ": no scenario 5-10 metrics yet - appendix figures skipped")
    next
  }
  bias_rmse_figure(o[["outcome"]], o[["prefix"]], scenarios = SUPP_SCENARIOS,
                   labels = SS_SUPP_SCENARIO_PLOT_LABELS, suffix = "_supp",
                   width = SUPP_WIDTH)
  rejection_figure(o[["outcome"]], o[["prefix"]], scenarios = SUPP_SCENARIOS,
                   labels = SS_SUPP_SCENARIO_PLOT_LABELS, suffix = "_supp",
                   width = SUPP_WIDTH)
}
