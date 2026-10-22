##########
# title: figures for the thesis appendix - crossfitting comparison
##########
# The appendix's case for the crossfitting strategy each method uses
# (crossfitting/), one panel per estimator family, the scenarios on x:
#   - cf_appendix_rmse_ratio.png: each arm's training-sample RMSE over the
#     double-crossfit arm's, per run (all arms of a run share one fold
#     assignment, so the ratio is paired), mean +/- 95% MCSE, line at 1;
#   - cf_appendix_ci.png: the CI pilot (crossfitting/confidence_intervals/),
#     half-sample bootstrap only - the one interval every arm has: marginal
#     coverage (line at 0.95) above mean interval length.
# Colour is the procedure, shared across families (CF_ARMS, R/figures.R), and
# the filled point is the arm the main analyses use. The table is
# thesis_tables/cf_tables.R.
#
# Writes to ../results/thesis_figures/crossfitting/. A study with no metrics
# file is skipped with a message.
#
#   cd results_processing/thesis_figures; Rscript crossfitting.R
#   Rscript crossfitting.R [results root] [figure dir]   # mock smoke run

library(here)
library(patchwork)
source(here("R", "figures.R"))

# paths
args <- commandArgs(trailingOnly = TRUE)
res_root <- if (length(args) >= 1) args[1] else file.path(dirname(here()), "results")
fig_path <- if (length(args) >= 2) {
  args[2]
} else {
  file.path(dirname(here()), "results", "thesis_figures", "crossfitting")
}
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

USED <- "Used in main analyses"
COMPARATOR <- "Comparator"

# named, so a procedure keeps its colour in a figure without the Double CF level
PROCEDURES <- unique(CF_ARMS$procedure)
PROCEDURE_COLOURS <- setNames(
  as.character(paletteer::paletteer_d("rcartocolor::Safe", length(PROCEDURES))),
  PROCEDURES
)

#' A metrics file under res_root, or NULL with a message if it isn't there
read_metrics <- function(rel_path) {
  f <- file.path(res_root, rel_path)
  if (!file.exists(f)) {
    message(f, " not found - its figure is skipped")
    return(NULL)
  }
  readRDS(f)
}

#' Join CF_ARMS and turn scenario, family, procedure and production into
#' display factors
label_arms <- function(df) {
  df %>%
    inner_join(CF_ARMS, by = c("arm", "family")) %>%
    mutate(
      scenario = label_factor(scenario, CF_SCENARIO_LABELS),
      family = label_factor(family, CF_FAMILY_LABELS),
      procedure = factor(procedure, levels = PROCEDURES),
      production = factor(ifelse(production, USED, COMPARATOR),
                          levels = c(USED, COMPARATOR))
    )
}

#' point_range_plot() in this script's layout: scenario on x, procedure by
#' colour, the production arm filled, one panel per family
cf_plot <- function(summary, metric, y_lab, hline, facet_scales = "fixed") {
  point_range_plot(
    summary,
    metric,
    y_lab,
    x = "scenario",
    colour = "procedure",
    shape = "production",
    facet_rows = NULL,
    facet_cols = "family",
    facet_scales = facet_scales,
    palette = scale_colour_manual(values = PROCEDURE_COLOURS),
    shape_palette = scale_shape_manual(values = setNames(c(16, 1), c(USED, COMPARATOR))),
    hline = hline,
    dodge_width = 0.6,
    point_size = 2
  ) +
    # wrapped rather than rotated: "Single effects + interaction" is long
    scale_x_discrete(labels = function(x) {
      vapply(x, function(s) paste(strwrap(s, 10), collapse = "\n"), "")
    }) +
    labs(x = NULL, colour = "Procedure", shape = NULL) +
    guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2)) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
}

# ---- RMSE relative to double crossfitting -----------------------------------

metrics <- read_metrics(file.path("crossfitting", "cf_metrics.RDS"))

if (!is.null(metrics)) {
  train <- metrics %>%
    filter(set == "train") %>%
    label_arms()

  ref <- train %>%
    filter(reference) %>%
    select(scenario, family, run, rmse_ref = rmse)

  ratio_summary <- train %>%
    filter(!reference) %>%
    inner_join(ref, by = c("scenario", "family", "run")) %>%
    mutate(rmse_ratio = rmse / rmse_ref) %>%
    summarise_metrics(c("scenario", "family", "procedure", "production"),
                      cols = c(rmse_ratio = "rmse_ratio"))

  fig <- cf_plot(ratio_summary, "rmse_ratio",
                 "RMSE / RMSE of double crossfitting\n(training sample)",
                 hline = 1) +
    theme(legend.position = "bottom", legend.box = "vertical")
  save_fig("cf_appendix_rmse_ratio.png", fig_path, height = 12, plot = fig)
}

# ---- CI pilot ----------------------------------------------------------------

ci_metrics <- read_metrics(file.path("crossfitting_ci", "cf_ci_metrics.RDS"))

if (!is.null(ci_metrics)) {
  ci_summary <- ci_metrics %>%
    filter(ci_method == "half_boot") %>%
    label_arms() %>%
    summarise_metrics(c("scenario", "family", "procedure", "production"),
                      cols = c(marg_cov = "marginal_coverage",
                               ci_len = "mean_ci_length"))

  cov_row <- cf_plot(ci_summary, "marg_cov", "Marginal coverage", hline = 0.95) +
    theme(axis.text.x = element_blank())

  # scenario names on the top row's strips already
  len_row <- cf_plot(ci_summary, "ci_len", "Mean interval length", hline = NULL,
                     facet_scales = "free_y") +
    theme(strip.text.x = element_blank(), strip.background.x = element_blank())

  fig <- (cov_row / len_row) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.box = "vertical")
  save_fig("cf_appendix_ci.png", fig_path, height = 18, plot = fig)
}
