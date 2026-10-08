##########
# title: figures for the thesis chapter - missing covariates, MI intervals
##########
# The MI confidence-interval example (missing/ci_example/): how to build a CATE
# interval after multiple imputation, continuous outcome only. One cell of the
# missing-data design - n = 500, 30% missing in both covariate types, MAR,
# multiple imputation - so what varies is scenario x estimator x the strategy
# combine_mi_ci() (R/bootstrap_ci.R) uses to pool the per-imputation bootstrap
# intervals. A separate script from missing.R, as ss_ci.R is from
# sample_size.R: different metrics file, metrics and design.
#
# One figure, three rows sharing the x axis: marginal coverage, simultaneous
# coverage and mean interval length, each by pooling strategy on x, estimator
# by colour, units by shape, scenario across. Notes:
#   - every metric is scored on all units, the complete units (no amputed
#     covariate) and the incomplete units (interval_metrics_split() in
#     R/metrics.R). An incomplete unit's interval is aimed at tau at its
#     unamputed covariates, which the data cannot pin down. MI keeps every
#     unit, so every estimator has all three;
#   - simultaneous coverage is 0/1 per run, so it gets the binomial MCSE,
#     marginal coverage and length the general one (summarise_metrics());
#   - every panel has its own y scale, as the scenarios differ. facet_grid()
#     can only free y per row, so each panel is its own plot, laid out 3 x 4
#     with patchwork, which keeps the columns aligned. Coverage panels have the
#     nominal 0.95 dashed, which their scale always takes in;
#   - each estimator keeps the colour it has in the other thesis figures (Safe,
#     by its place in MODEL_LABELS, as ss_ci.R), although only four of them
#     (CI_MODELS) have intervals.
#
# The scenarios are the main study's 1-4 (missing/ci_example/README.md), so
# they take the sample-size chapter's labels, as missing.R does.
#
# Writes to ../results/thesis_figures/miss_cts_ci/miss_cts_ci_all.png

library(here)
library(patchwork)
source(here("R", "figures.R"))

# paths
path <- here()
res_path <- file.path(dirname(path), "results", "missing", "ci_example")
fig_path <- file.path(dirname(path), "results", "thesis_figures", "miss_cts_ci")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

NOMINAL <- 0.95

# one row of the figure each, top to bottom; `binomial` marks a per-run 0/1
# indicator, `coverage` gets the nominal line
CI_METRICS <- tibble::tribble(
  ~stem,       ~col,                    ~lab,                    ~coverage, ~binomial,
  "marg_cov",  "marginal_coverage",     "Marginal coverage",     TRUE,      FALSE,
  "simul_cov", "simultaneous_coverage", "Simultaneous coverage", TRUE,      TRUE,
  "ci_len",    "mean_ci_length",        "Mean interval length",  FALSE,     FALSE
)

# as missing.R
UNIT_LABELS <- c(all = "All participants", cu = "Complete participants",
                 iu = "Incomplete participants")
UNIT_SHAPES <- setNames(c(16, 17, 15), UNIT_LABELS)

metrics <- readRDS(file.path(res_path, "cts_miss_ci_metrics.RDS")) %>%
  filter(scenario %in% 1:4)

# one row per run and set of units: the all-unit columns have no suffix, so
# give them one before splitting the names
metrics_long <- metrics %>%
  select(scenario, model, strategy, run,
         all_of(c(CI_METRICS$col, paste0(CI_METRICS$col, "_cu"),
                  paste0(CI_METRICS$col, "_iu")))) %>%
  rename_with(~ paste0(.x, "_all"), all_of(CI_METRICS$col)) %>%
  pivot_longer(-c(scenario, model, strategy, run),
               names_to = c(".value", "units"),
               names_pattern = "^(.*)_(all|cu|iu)$") %>%
  mutate(units = factor(UNIT_LABELS[units], levels = UNIT_LABELS))

metrics_summary <- summarise_metrics(
  metrics_long,
  c("scenario", "model", "strategy", "units"),
  cols = setNames(CI_METRICS$col, CI_METRICS$stem),
  binomial = CI_METRICS$stem[CI_METRICS$binomial],
  count_na = character()
) %>%
  apply_labels(SS_SCENARIO_LABELS)

model_levels <- levels(metrics_summary$model)
safe <- as.character(paletteer_d("rcartocolor::Safe", length(MODEL_LABELS)))
pal_values <- setNames(safe[match(model_levels, MODEL_LABELS)], model_levels)

#' One panel: one metric in one scenario
#'
#' Up to 12 points per strategy (4 estimators x 3 sets of units), so the dodge
#' is wide, the points bigger and the error bars fainter, as missing.R's
#' models_figure().
#'
#' @param m one row of CI_METRICS
#' @param scen a scenario label
#' @param top TRUE for the top row, the only one with the scenario strips
#' @param bottom TRUE for the bottom row, the only one with the x axis labels
#'   and title
#' @param left TRUE for the first column, the only one with the y title
ci_panel <- function(m, scen, top = FALSE, bottom = FALSE, left = FALSE) {
  p <- point_range_plot(
    filter(metrics_summary, scenario == scen),
    m$stem,
    m$lab,
    x = "strategy",
    colour = "model",
    shape = "units",
    shape_palette = scale_shape_manual(values = UNIT_SHAPES),
    facet_rows = NULL,
    facet_cols = "scenario",
    hline = if (m$coverage) NOMINAL else NULL,
    palette = scale_colour_manual(values = pal_values, limits = model_levels),
    dodge_width = 0.85,
    point_size = 1.5,
    ci_alpha = 0.5
  ) +
    labs(x = "Pooling strategy", colour = NULL, shape = NULL) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  if (!top) {
    p <- p + theme(strip.text.x = element_blank(),
                   strip.background.x = element_blank())
  }
  if (!bottom) {
    p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  }
  if (!left) {
    p <- p + theme(axis.title.y = element_blank())
  }
  p
}

scenarios <- levels(droplevels(metrics_summary$scenario))
n_rows <- nrow(CI_METRICS)
panels <- list()
for (i in seq_len(n_rows)) {
  for (j in seq_along(scenarios)) {
    panels[[length(panels) + 1]] <- ci_panel(
      CI_METRICS[i, ], scenarios[j],
      top = i == 1, bottom = i == n_rows, left = j == 1
    )
  }
}

fig <- wrap_plots(panels, ncol = length(scenarios)) +
  plot_layout(guides = "collect", axis_titles = "collect_x") &
  theme(legend.position = "bottom", legend.box = "vertical")

save_fig("miss_cts_ci_all.png", fig_path, width = 30, height = 24, plot = fig)
