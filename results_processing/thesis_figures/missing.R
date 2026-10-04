##########
# title: figures for the thesis chapter - missing covariates, continuous outcome
##########
# The main chapter's missing-data figures, continuous outcome only, for one
# estimator: the DR-learner with random forests. Both metrics compare each
# handling method with the complete-data arm (no amputation), same scenario,
# mechanism and run:
#   - bias difference, bias - complete-data bias, per run and averaged over
#     runs: the bias the missingness and its handling add. 0 is no change.
#     The same measure as bias_diff_complete_cu in cts_miss_metrics.R and
#     missing/continuous/cts_miss_results.qmd, extended to the incomplete units;
#   - relative efficiency, MSE / complete-data MSE (rel_efficiency_cu/_iu from
#     cts_miss_metrics.R, a per-run ratio averaged over runs, as in the
#     appendix tables). 1 is no loss.
# One figure each, the mechanism x scenario grid with handling method on x.
# The complete-data arm is 0 and 1 by construction, so it is not plotted.
# (All five estimators at once was too busy - the appendix tables,
# thesis_tables/miss_tables.R, have the rest.)
#
# Every metric is scored separately on the complete units (no amputed
# covariate) and the incomplete units (cate_metrics_split() in R/metrics.R),
# shown in the same panel by colour. The complete-data arm is split by the
# amputation its run would have had, so each set of units is compared with the
# same units in the reference. Notes:
#   - the truth is the primary one, tau(X), under MNAR-tau too;
#   - complete case and IPW drop the incomplete units, so they have only a
#     complete-unit point;
#   - scenario 1 has no MNAR-tau arm (cts_miss_config.R), so that panel is
#     empty.
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
REF_METHOD <- "complete_data"
UNIT_LABELS <- c(cu = "Complete units", iu = "Incomplete units")

metrics <- readRDS(file.path(res_path, "cts_miss_metrics.RDS"))

# one row per run and set of units
metrics_long <- metrics %>%
  filter(scenario %in% 1:4, model == MODEL) %>%
  select(scenario, mechanism, method, run,
         bias_cu, bias_iu, rel_efficiency_cu, rel_efficiency_iu) %>%
  pivot_longer(c(bias_cu, bias_iu, rel_efficiency_cu, rel_efficiency_iu),
               names_to = c(".value", "units"),
               names_pattern = "^(bias|rel_efficiency)_(cu|iu)$")

# the complete-data arm's bias on the same run and units. cts_miss_metrics.R
# keeps the difference for the complete units only (bias_diff_complete_cu), so
# take it for both sets of units here; on the complete units it is the same
# number. Complete case and IPW's incomplete-unit bias is all NA, so theirs
# comes out NA too.
ref_bias <- metrics_long %>%
  filter(method == REF_METHOD) %>%
  select(scenario, mechanism, run, units, bias_ref = bias)

metrics_long <- metrics_long %>%
  filter(method != REF_METHOD) %>%
  left_join(ref_bias, by = c("scenario", "mechanism", "run", "units")) %>%
  mutate(bias_diff = bias - bias_ref) %>%
  mutate(units = factor(UNIT_LABELS[units], levels = UNIT_LABELS)) %>%
  apply_labels(SS_SCENARIO_LABELS)

metrics_summary <- summarise_metrics(
  metrics_long,
  c("scenario", "mechanism", "method", "units"),
  cols = c(bias_diff = "bias_diff", rel_efficiency = "rel_efficiency"),
  count_na = character()
)

#' One metric by handling method, complete and incomplete units by colour
#'
#' @param summary metrics_summary
#' @param metric "bias_diff" or "rel_efficiency"
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
    theme(legend.position = "bottom") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  save_fig(file, fig_path, height = 18, plot = fig)
  fig
}

bias_fig <- units_figure(metrics_summary, "bias_diff",
                         "Bias minus complete-data bias",
                         "cts_miss_bias_diff.png")
eff_fig <- units_figure(metrics_summary, "rel_efficiency",
                        "Relative efficiency vs complete data",
                        "cts_miss_rel_efficiency.png", hline = 1)
