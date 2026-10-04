##########
# title: figures for the thesis chapter - missing covariates, continuous outcome
##########
# The main chapter's missing-data figures, continuous outcome only, for one
# estimator: the DR-learner with random forests. CATE bias in one figure, RMSE
# in another. Each is the mechanism x scenario grid, with handling method on x.
# (All five estimators at once was too busy - the appendix tables,
# thesis_tables/miss_tables.R, have the rest.)
#
# Every metric is scored separately on the complete units (no amputed
# covariate) and the incomplete units (cate_metrics_split() in R/metrics.R),
# shown in the same panel by colour. The incomplete units carry the error
# floor of scoring against tau at covariates no method saw (missing/ADEMP.md,
# "Estimands"). Notes:
#   - the truth is the primary one, tau(X), under MNAR-tau too;
#   - complete case and IPW drop the incomplete units, so they have only a
#     complete-unit point;
#   - scenario 1 has no MNAR-tau arm (cts_miss_config.R), so that panel is
#     empty;
#   - the complete-data arm is the reference, plotted as a method of its own,
#     scored on the units its run's amputation would have made incomplete.
#
# The main scenarios only, 1-4. These are the main study's scenarios 1-4
# (TE_MISS in R/dgm_scenarios.R), so they take the sample-size chapter's labels.
# Labels, palette, summaries and figure sizing come from R/figures.R. This
# script carries only the paths and this study's filters.

library(here)
source(here("R", "figures.R"))

# paths
path <- here()
res_path <- file.path(dirname(path), "results", "missing", "continuous")
fig_path <- file.path(dirname(path), "results", "thesis_figures", "missing")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

MODEL <- "dr_random_forest"
UNIT_LABELS <- c(cu = "Complete units", iu = "Incomplete units")

metrics <- readRDS(file.path(res_path, "cts_miss_metrics.RDS"))

# one row per run and set of units
metrics_long <- metrics %>%
  filter(scenario %in% 1:4, model == MODEL) %>%
  select(scenario, mechanism, method, run,
         bias_cu, bias_iu, rmse_cu, rmse_iu) %>%
  pivot_longer(c(bias_cu, bias_iu, rmse_cu, rmse_iu),
               names_to = c(".value", "units"),
               names_pattern = "^(bias|rmse)_(cu|iu)$") %>%
  mutate(units = factor(UNIT_LABELS[units], levels = UNIT_LABELS)) %>%
  apply_labels(SS_SCENARIO_LABELS)

metrics_summary <- summarise_metrics(
  metrics_long,
  c("scenario", "mechanism", "method", "units"),
  cols = c(bias = "bias", rmse = "rmse"),
  count_na = character()
)

#' One metric by handling method, complete and incomplete units by colour
#'
#' @param summary metrics_summary
#' @param metric "bias" or "rmse"
#' @param y_lab axis label
#' @param file output file name
#' @param hline reference line, as point_range_plot()
units_figure <- function(summary, metric, y_lab, file, hline = 0) {
  # complete case and IPW's incomplete-unit rows are all NA (mean NaN), so
  # leave them out rather than give them a dodge slot
  keep <- filter(summary, is.finite(.data[[paste0("mean_", metric)]]))

  fig <- point_range_plot(keep, metric, y_lab, x = "method", colour = "units",
                          hline = hline) +
    labs(x = "Missing-data handling", colour = NULL) +
    theme(legend.position = "bottom")
  save_fig(file, fig_path, height = 18, plot = fig)
  fig
}

bias_fig <- units_figure(metrics_summary, "bias",
                         "Bias of the CATE (DR-RandomForest)",
                         "cts_miss_bias.png")
rmse_fig <- units_figure(metrics_summary, "rmse",
                         "RMSE of the CATE (DR-RandomForest)",
                         "cts_miss_rmse.png", hline = NULL)
