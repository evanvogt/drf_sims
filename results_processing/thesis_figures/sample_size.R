##########
# title: figures for the thesis chapter - sample size, correlated covariates
##########
# The main chapter's sample-size figures: the correlated-covariate studies
# (sample_size/correlated/) at rho = 0.5. The independent-covariate studies
# (cts_ss.R, bin_ss.R) are the appendix's.
#
# One figure per outcome: CATE bias on the top row, RMSE on the bottom, the
# four scenarios across. Labels, palette, summaries and figure sizing come from
# R/figures.R. This script carries only the paths and this study's filters.

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
