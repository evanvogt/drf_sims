##########
# title: figures for the thesis chapter - missing covariates
##########
# The main chapter's missing-data figures, continuous and binary outcomes. Both
# metrics compare each handling method with the complete-data arm (no
# amputation), same scenario, mechanism, estimator and run:
#   - bias difference, bias - complete-data bias, per run and averaged over
#     runs: the bias the missingness and its handling add. 0 is no change.
#     The same measure as bias_diff_complete/_cu in missing/*/*_miss_metrics.R
#     and missing/*/*_miss_results.qmd, extended to the incomplete units;
#   - relative efficiency, MSE / complete-data MSE (rel_efficiency/_cu/_iu from
#     *_miss_metrics.R, a per-run ratio averaged over runs, as in the
#     appendix tables). 1 is no loss.
# Each is drawn on the mechanism x scenario grid with handling method on x, for
# each outcome (prefix cts / bin):
#   - one figure per estimator, units by colour
#     (<prefix>_miss_<metric>_<model>.png);
#   - one figure with every estimator, estimator by colour and units by shape
#     (<prefix>_miss_<metric>_all_models.png).
# The complete-data arm is 0 and 1 by construction, so it is not plotted. The
# binary outcome's CATE is a risk difference, so its bias difference is on that
# scale.
#
# Every metric is scored on all units, and separately on the complete units (no
# amputed covariate) and the incomplete units (cate_metrics_split() in
# R/metrics.R) - the same three columns as the appendix tables. The
# complete-data arm is split by the amputation its run would have had, so each
# set of units is compared with the same units in the reference. Notes:
#   - the truth is the primary one, tau(X), under MNAR-tau too;
#   - complete case and IPW drop the incomplete units, so they have only a
#     complete-unit point (*_miss_metrics.R leaves their all-unit metrics NA:
#     a ~350-unit MSE over a 500-unit one is not a like-for-like ratio);
#   - not every estimator has every handling method in the metrics files (in
#     the continuous results MIA is causal forest and DR-RandomForest only, and
#     MI is not DR-oracle or DR-SuperLearner), so those slots are empty;
#   - scenario 1 has no MNAR-tau arm (missing/*/*_miss_config.R), so that panel
#     is empty.
#
# The main scenarios only, 1-4. These are the main study's scenarios 1-4
# (TE_MISS in R/dgm_scenarios.R), so they take the sample-size chapter's labels.
# Labels, palette, summaries and figure sizing come from R/figures.R. This
# script carries only the paths and this study's filters.

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
