##########
# title: figures for the thesis chapter - missing covariates
##########
# Main-chapter missing-data figures, continuous and binary, scenarios 1-4. Each
# handling method is compared with the complete-data arm (same scenario,
# mechanism, estimator and run): bias difference (0 = no change) and relative
# efficiency, MSE / complete-data MSE (1 = no loss), on all, complete and
# incomplete units, as in the appendix tables.
#
# - Complete case and IPW have only a complete-unit point (their all-unit ratio
#   would compare ~350 units with 500).
# - Empty slots are real: not every estimator has every handling method, and
#   scenario 1 has no MNAR-tau arm.
# - The truth is tau(X), under MNAR-tau too; the binary bias is on the
#   risk-difference scale.
#
# Writes, per outcome prefix (cts, bin): <prefix>_miss_<metric>_<model>.png
# (units by colour) and <prefix>_miss_<metric>_all_models.png.

library(here)
source(here("R", "figures.R"))

# paths
path <- here()
res_root <- file.path(dirname(path), "results", "missing")
fig_path <- file.path(dirname(path), "results", "thesis_figures", "missing")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

outcomes <- list(
  list(dir = "continuous", prefix = "cts"),
  list(dir = "binary", prefix = "bin")
)

REF_METHOD <- "complete_data"
UNIT_LABELS <- c(all = "All participants", cu = "Complete participants",
                 iu = "Incomplete participants")
UNIT_SHAPES <- setNames(c(16, 17, 15), UNIT_LABELS)

FIGURES <- list(
  bias_diff = list(y_lab = "Bias minus complete-data bias", hline = 0),
  rel_efficiency = list(y_lab = "Relative efficiency vs complete data",
                        hline = 1)
)

#' Mean and MCSE of both metrics, by estimator, method and set of units
#'
#' @param metrics one outcome's *_miss_metrics.RDS
summarise_units <- function(metrics) {
  metrics <- filter(metrics, scenario %in% 1:4)

  # the metrics file has the bias difference on all units and the complete
  # units only, so take the incomplete units' here, as
  # thesis_tables/miss_tables.R does. Complete case and IPW's incomplete-unit
  # bias is all NA, so theirs comes out NA too.
  ref_bias_iu <- metrics %>%
    filter(method == REF_METHOD) %>%
    select(scenario, mechanism, model, run, bias_iu_complete = bias_iu)

  # one row per run and set of units
  metrics_long <- metrics %>%
    filter(method != REF_METHOD) %>%
    left_join(ref_bias_iu, by = c("scenario", "mechanism", "model", "run")) %>%
    mutate(bias_diff_iu = bias_iu - bias_iu_complete) %>%
    select(scenario, mechanism, model, method, run,
           bias_diff_all = bias_diff_complete,
           bias_diff_cu = bias_diff_complete_cu, bias_diff_iu,
           rel_efficiency_all = rel_efficiency,
           rel_efficiency_cu, rel_efficiency_iu) %>%
    pivot_longer(-c(scenario, mechanism, model, method, run),
                 names_to = c(".value", "units"),
                 names_pattern = "^(bias_diff|rel_efficiency)_(all|cu|iu)$") %>%
    mutate(units = factor(UNIT_LABELS[units], levels = UNIT_LABELS)) %>%
    apply_labels(SS_SCENARIO_LABELS)

  # complete case and IPW's all- and incomplete-unit rows are all NA (mean
  # NaN), as is every estimator x method pair that wasn't run, so the figures
  # leave out non-finite means rather than give them a dodge slot
  summarise_metrics(
    metrics_long,
    c("scenario", "mechanism", "model", "method", "units"),
    cols = c(bias_diff = "bias_diff", rel_efficiency = "rel_efficiency"),
    count_na = character()
  )
}

#' One metric for one estimator, all, complete and incomplete units by colour
#'
#' @param summary summarise_units() output, filtered to one model
#' @param metric "bias_diff" or "rel_efficiency"
#' @param y_lab axis label
#' @param file output file name
#' @param hline reference line, as point_range_plot()
units_figure <- function(summary, metric, y_lab, file, hline = 0) {
  keep <- filter(summary, is.finite(.data[[paste0("mean_", metric)]]))

  fig <- point_range_plot(keep, metric, y_lab, x = "method", colour = "units",
                          hline = hline) +
    labs(x = "Missing-data handling", colour = NULL) +
    theme(legend.position = "bottom") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  save_fig(file, fig_path, height = 18, plot = fig)
  fig
}

#' One metric for every estimator: estimator by colour, units by shape
#'
#' Up to 15 points per handling method (5 estimators x 3 sets of units), so the
#' dodge is wider, the points bigger and the error bars fainter than
#' units_figure(), and the figure is bigger.
#'
#' @inheritParams units_figure
#' @param summary summarise_units() output, all models
models_figure <- function(summary, metric, y_lab, file, hline = 0) {
  keep <- filter(summary, is.finite(.data[[paste0("mean_", metric)]]))

  fig <- point_range_plot(keep, metric, y_lab, x = "method", colour = "model",
                          shape = "units",
                          shape_palette = scale_shape_manual(values = UNIT_SHAPES),
                          hline = hline, dodge_width = 0.85, point_size = 1.5,
                          ci_alpha = 0.5) +
    labs(x = "Missing-data handling", colour = NULL, shape = NULL) +
    theme(legend.position = "bottom", legend.box = "vertical") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  save_fig(file, fig_path, width = 30, height = 24, plot = fig)
  fig
}

#' Every figure for one outcome
#'
#' @param o an element of `outcomes`
outcome_figures <- function(o) {
  metrics <- readRDS(file.path(res_root, o$dir,
                               paste0(o$prefix, "_miss_metrics.RDS")))
  metrics_summary <- summarise_units(metrics)

  # one figure per estimator and metric, named by the estimator's raw name
  models <- intersect(names(MODEL_LABELS), unique(metrics$model))
  model_figs <- lapply(setNames(models, models), function(m) {
    summary <- filter(metrics_summary, model == MODEL_LABELS[[m]])
    Map(function(metric, spec) {
      units_figure(summary, metric, spec$y_lab,
                   paste0(o$prefix, "_miss_", metric, "_", m, ".png"),
                   hline = spec$hline)
    }, names(FIGURES), FIGURES)
  })

  all_models <- Map(function(metric, spec) {
    models_figure(metrics_summary, metric, spec$y_lab,
                  paste0(o$prefix, "_miss_", metric, "_all_models.png"),
                  hline = spec$hline)
  }, names(FIGURES), FIGURES)

  list(by_model = model_figs, all_models = all_models)
}

figs <- lapply(setNames(outcomes, sapply(outcomes, `[[`, "prefix")),
               outcome_figures)
