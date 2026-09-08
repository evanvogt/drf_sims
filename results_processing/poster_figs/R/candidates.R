##########
# title: candidate plot builders for the ICTMC poster - pilot panels
##########
# Exploration only, not the final poster figures: for the sample-size and
# missing-data panels, several metric/encoding choices are still open (see
# ../ICTMC_poster_planning.md). These functions build one lollipop_plot() per
# candidate encoding so explore_gallery.R and explore_options.qmd can show
# them all side by side without duplicating the plot-building logic.
#
# Piloted on these two panels only, because they already have working data
# prep to build from (results_processing/ictmc_figs.R and
# results_processing/thesis_figures/{cts_ss,bin_ss,miss_cts,miss_bin}.R). The
# other four poster panels (HTE testing, CI coverage, model evaluation,
# subgroup validation) still need that data prep written before this pattern
# can extend to them.

library(here)
library(dplyr)
source(here("R", "figures.R"))

# ---- data ---------------------------------------------------------------

metrics_all <- readRDS(file.path(
  dirname(here()),
  "collected_metrics",
  "all_metrics.RDS"
))

# same filters as ictmc_figs.R's metrics_ss, but summarised via
# summarise_metrics()/apply_labels() instead of a hand-rolled group_by block
metrics_ss_raw <- metrics_all %>%
  filter(study_name %in% c("binary", "continuous")) %>%
  filter(model != "dr_oracle") %>%
  filter(scenario %in% c(1, 3, 8, 9)) %>%
  apply_labels(SS_SCENARIO_LABELS) %>%
  mutate(n = factor(n, levels = c(100, 250, 500, 1000))) %>%
  droplevels()

# same filters as ictmc_figs.R's metrics_miss - no scenario subset, keeps
# whichever scenarios the missing-data studies ran (MISS_SCENARIO_LABELS'
# five-scenario set)
metrics_miss_raw <- metrics_all %>%
  filter(study_name %in% c("missing/binary", "missing/continuous")) %>%
  filter(model != "dr_oracle") %>%
  apply_labels(MISS_SCENARIO_LABELS) %>%
  droplevels()

# ---- summaries ------------------------------------------------------------

# rmse and spearman aren't in summarise_metrics()'s default `cols` - both are
# candidate metrics for this panel per the planning doc
SS_METRIC_COLS <- c(
  bias = "bias", rel_ate_bias = "rel_ate_bias", rel_bias_cate = "rel_bias_cate",
  mse = "mse", rmse = "rmse", spearman = "spearman"
)

# rel_efficiency and rel_bias_complete (missing/*/*_miss_metrics.R) are ratios
# against the complete-data arm - see R/metrics.R's cate_metrics() roxygen for
# how these differ from rel_ate_bias/rel_bias_cate
MISS_METRIC_COLS <- c(
  bias = "bias", rel_bias_cate = "rel_bias_cate", mse = "mse", rmse = "rmse",
  rel_efficiency = "rel_efficiency", rel_bias_complete = "rel_bias_complete"
)

# kept separate by study_name (outcome type) so candidates can facet on it
# instead of averaging over it - see candidate_plots_sample_size()'s `outcome`
ss_summary <- summarise_metrics(
  metrics_ss_raw,
  c("scenario", "n", "model", "study_name"),
  cols = SS_METRIC_COLS
)

# averaged over CATE model - the summarisation the planning doc leans towards,
# to keep the missing-data plot to one line per handling method
miss_summary_avg <- summarise_metrics(
  metrics_miss_raw,
  c("scenario", "mechanism", "method", "study_name"),
  cols = MISS_METRIC_COLS
)

# kept split by model, to check whether averaging hides a method x model
# interaction worth keeping
miss_summary_by_model <- summarise_metrics(
  metrics_miss_raw,
  c("scenario", "mechanism", "method", "model", "study_name"),
  cols = MISS_METRIC_COLS
)

# ---- candidate builders -----------------------------------------------------

#' Candidate encodings for the sample-size panel
#'
#' One lollipop per metric raised as an option in the planning doc: bias,
#' relative ATE bias, relative CATE bias, MSE, RMSE, and Spearman correlation
#' (shown as deviation from 1, not 0, since lower correlation is worse). x is
#' sample size, colour is model, faceted by scenario - the shape the
#' sample-size studies actually use (thesis_figures/cts_ss.R), not
#' point_range_plot()'s default mechanism x scenario grid.
#'
#' @param outcome "both" facets outcome type as an extra row (study_name);
#'   "binary"/"continuous" restrict to one outcome and drop that facet - the
#'   three-way choice the planning doc raises for this panel
candidate_plots_sample_size <- function(summary = ss_summary,
                                        outcome = c("both", "binary", "continuous")) {
  outcome <- match.arg(outcome)
  if (outcome != "both") summary <- filter(summary, study_name == outcome)

  common <- list(
    x = "n", colour = "model", facet_cols = "scenario",
    facet_rows = if (outcome == "both") "study_name" else NULL
  )

  specs <- list(
    bias          = list(metric = "bias",          y_lab = "Mean bias",                    hline = 0),
    rel_ate_bias  = list(metric = "rel_ate_bias",   y_lab = "Mean relative ATE bias",        hline = 0),
    rel_bias_cate = list(metric = "rel_bias_cate",  y_lab = "Mean relative CATE bias",       hline = 0),
    mse           = list(metric = "mse",            y_lab = "Mean MSE",                      hline = 0),
    rmse          = list(metric = "rmse",           y_lab = "Mean RMSE",                     hline = 0),
    spearman_dev  = list(metric = "spearman",       y_lab = "Mean Spearman correlation",     hline = 1)
  )

  lapply(specs, function(s) {
    do.call(lollipop_plot, c(list(summary, s$metric, s$y_lab, hline = s$hline), common))
  })
}

#' Candidate encodings for the missing-data panel
#'
#' One lollipop per metric raised as an option in the planning doc: bias,
#' relative CATE bias, MSE, RMSE, relative efficiency and relative-to-complete
#' bias (the latter two are ratios against the complete-data arm, so their
#' reference line is 1, not 0). x is missing-data mechanism, colour is
#' handling method, faceted by scenario.
#'
#' @param by_model FALSE (default) averages over CATE model, one line per
#'   method - the planning doc's preferred summarisation. TRUE keeps model as
#'   a further row facet, to compare against the averaged version.
candidate_plots_missing <- function(by_model = FALSE) {
  summary <- if (by_model) miss_summary_by_model else miss_summary_avg

  common <- list(
    x = "mechanism", colour = "method", facet_cols = "scenario",
    facet_rows = if (by_model) "model" else "study_name"
  )

  specs <- list(
    bias              = list(metric = "bias",              y_lab = "Mean bias",                                        hline = 0),
    rel_bias_cate     = list(metric = "rel_bias_cate",      y_lab = "Mean relative CATE bias",                          hline = 0),
    mse               = list(metric = "mse",                y_lab = "Mean MSE",                                         hline = 0),
    rmse              = list(metric = "rmse",               y_lab = "Mean RMSE",                                        hline = 0),
    rel_efficiency    = list(metric = "rel_efficiency",     y_lab = "Relative efficiency (MSE / complete-data MSE)",    hline = 1),
    rel_bias_complete = list(metric = "rel_bias_complete",  y_lab = "Relative bias (bias / complete-data bias)",        hline = 1)
  )

  lapply(specs, function(s) {
    do.call(lollipop_plot, c(list(summary, s$metric, s$y_lab, hline = s$hline), common))
  })
}
