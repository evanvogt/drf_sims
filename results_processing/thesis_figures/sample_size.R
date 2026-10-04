##########
# title: figures for the thesis chapter - sample size, correlated covariates
##########
# The main chapter's sample-size figures: the correlated-covariate studies
# (sample_size/correlated/) at rho = 0.5. The independent-covariate studies
# (cts_ss.R, bin_ss.R) are the appendix's.
#
# Two figures per outcome, the four scenarios across in each:
#   - CATE bias on the top row, RMSE on the bottom;
#   - each HTE test's rejection rate at 0.05, one row per test, with the same
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

#' Bias (top) and RMSE (bottom) by sample size, one panel per scenario
#'
#' @param outcome results folder, "continuous" or "binary"
#' @param prefix file prefix, "cts" or "bin"
#' @param scale_lab appended to the y labels, e.g. " (risk difference)"
bias_rmse_figure <- function(outcome, prefix, scale_lab = "") {
  metrics <- readRDS(file.path(res_path, outcome,
                               paste0(prefix, "_corr_metrics.RDS"))) %>%
    filter(rho == RHO) %>%
    apply_labels(SS_SCENARIO_LABELS) %>%
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
  save_fig(paste0(prefix, "_corr_bias_rmse.png"), fig_path, plot = fig)
  fig
}

cts_fig <- bias_rmse_figure("continuous", "cts")
bin_fig <- bias_rmse_figure("binary", "bin")

# the HTE tests, in facet-row order
# the HTE tests, in facet-row order
HTE_TESTS <- c(
  BLP_p_os = "BLP (one-sided, HC3)",
  BLP_p = "BLP (two-sided)",
  indep_po = "PO independence",
  indep_cate = "CATE independence"
)

TRUE_LAB <- "True values"

#' Rejection rate at 0.05 by sample size, one row per test, one column per
#' scenario
#'
#' The "True values" series is the BLP and CATE permutation tests run on the
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
rejection_figure <- function(outcome, prefix) {
  read_corr <- function(stem) {
    readRDS(file.path(res_path, outcome, paste0(prefix, "_corr_", stem, ".RDS"))) %>%
      filter(rho == RHO)
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
    mutate(rej = as.numeric(p < 0.1),
           test = factor(test, levels = names(HTE_TESTS), labels = HTE_TESTS)) %>%
    apply_labels(SS_SCENARIO_LABELS) %>%
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

  fig <- point_range_plot(rej_summary, "rej", "Rejection rate at 0.1",
                          x = "n", colour = "model", facet_rows = "test",
                          facet_cols = "scenario", facet_scales = "fixed",
                          line = TRUE, ci_alpha = 0.7, hline = 0.1,
                          palette = pal) +
    labs(x = "Sample size", colour = "Model") +
    theme(legend.position = "right") #+ # figures are too long with legend at the bottom
    #guides(colour = guide_legend(nrow = 2))
  save_fig(paste0(prefix, "_corr_rejection.png"), fig_path, height = 18, plot = fig)
  fig
}

cts_rej_fig <- rejection_figure("continuous", "cts")
bin_rej_fig <- rejection_figure("binary", "bin")
