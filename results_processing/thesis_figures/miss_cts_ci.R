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
# by colour, scenario across. Notes:
#   - all units only. The metrics file also scores the complete and incomplete
#     units (_cu / _iu, interval_metrics_split() in R/metrics.R); the
#     diagnostic figures (missing/ci_example/cts_miss_ci_results.R) show those;
#   - simultaneous coverage is 0/1 per run, so it gets the binomial MCSE,
#     marginal coverage and length the general one (summarise_metrics());
#   - coverage is drawn on [0, 1] with the nominal 0.95 dashed; length shares
#     one y scale across the scenarios;
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
# indicator, `coverage` gets the nominal line and the [0, 1] scale
CI_METRICS <- tibble::tribble(
  ~stem,       ~col,                    ~lab,                    ~coverage, ~binomial,
  "marg_cov",  "marginal_coverage",     "Marginal coverage",     TRUE,      FALSE,
  "simul_cov", "simultaneous_coverage", "Simultaneous coverage", TRUE,      TRUE,
  "ci_len",    "mean_ci_length",        "Mean interval length",  FALSE,     FALSE
)

metrics <- readRDS(file.path(res_path, "cts_miss_ci_metrics.RDS")) %>%
  filter(scenario %in% 1:4)

metrics_summary <- summarise_metrics(
  metrics,
  c("scenario", "model", "strategy"),
  cols = setNames(CI_METRICS$col, CI_METRICS$stem),
  binomial = CI_METRICS$stem[CI_METRICS$binomial],
  count_na = character()
) %>%
  apply_labels(SS_SCENARIO_LABELS)

model_levels <- levels(metrics_summary$model)
safe <- as.character(paletteer_d("rcartocolor::Safe", length(MODEL_LABELS)))
pal_values <- setNames(safe[match(model_levels, MODEL_LABELS)], model_levels)

#' One row of the figure: one metric, scenarios across
#'
#' @param m one row of CI_METRICS
#' @param top TRUE for the top row, the only one with the scenario strips
#' @param bottom TRUE for the bottom row, the only one with the x axis labels
#'   and title
ci_row <- function(m, top = FALSE, bottom = FALSE) {
  p <- point_range_plot(
    metrics_summary,
    m$stem,
    m$lab,
    x = "strategy",
    colour = "model",
    facet_rows = NULL,
    facet_cols = "scenario",
    facet_scales = "fixed",
    hline = if (m$coverage) NOMINAL else NULL,
    palette = scale_colour_manual(values = pal_values, limits = model_levels)
  ) +
    labs(x = "Pooling strategy", colour = NULL) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  # coord, not scale, limits: a scale limit would drop an error bar that
  # crosses 0 or 1 rather than clip it
  if (m$coverage) p <- p + coord_cartesian(ylim = c(0, 1))
  if (!bottom) {
    p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  }
  if (!top) {
    p <- p + theme(strip.text.x = element_blank(),
                   strip.background.x = element_blank())
  }
  p
}

n_rows <- nrow(CI_METRICS)
rows <- lapply(seq_len(n_rows), function(i) {
  ci_row(CI_METRICS[i, ], top = i == 1, bottom = i == n_rows)
})

fig <- wrap_plots(rows, ncol = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

save_fig("miss_cts_ci_all.png", fig_path, width = 21, height = 20, plot = fig)
