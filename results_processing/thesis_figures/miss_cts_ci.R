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
#   - every panel has its own y scale, as the scenarios differ
#     (ggh4x::facet_grid2(independent = "y"), as facet_grid() can only free y
#     per row). The metric names are the row strips, on the left where a y
#     title would be. Coverage panels have the nominal 0.95 dashed, which
#     their scale always takes in;
#   - each estimator keeps the colour it has in the other thesis figures (Safe,
#     by its place in MODEL_LABELS, as ss_ci.R), although only four of them
#     (CI_MODELS) have intervals.
#
# The scenarios are the main study's 1-4 (missing/ci_example/README.md), so
# they take the sample-size chapter's labels, as missing.R does.
#
# Writes to ../results/thesis_figures/miss_cts_ci/miss_cts_ci_all.png

library(here)
library(ggh4x)
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

# one row per run of the figure: the metric is a facet row, its mean and MCSE
# in mean_value / mcse_value for point_range_plot()
metrics_rows <- metrics_summary %>%
  pivot_longer(matches("^(mean|mcse)_"), names_to = c(".value", "stem"),
               names_pattern = "^(mean|mcse)_(.*)$") %>%
  rename(mean_value = mean, mcse_value = mcse) %>%
  mutate(metric = factor(CI_METRICS$lab[match(stem, CI_METRICS$stem)],
                         levels = CI_METRICS$lab))

# the nominal line, in the coverage rows only
nominal <- tibble(
  metric = factor(CI_METRICS$lab[CI_METRICS$coverage], levels = CI_METRICS$lab),
  yintercept = NOMINAL
)

# Up to 12 points per strategy (4 estimators x 3 sets of units), so the dodge
# is wide, the points bigger and the error bars fainter, as missing.R's
# models_figure(). facet_grid2() replaces point_range_plot()'s facet_grid().
fig <- point_range_plot(
  metrics_rows,
  "value",
  NULL,
  x = "strategy",
  colour = "model",
  shape = "units",
  shape_palette = scale_shape_manual(values = UNIT_SHAPES),
  facet_rows = "metric",
  hline = nominal,
  palette = scale_colour_manual(values = pal_values, limits = model_levels),
  dodge_width = 0.85,
  point_size = 1.5,
  ci_alpha = 0.5
) +
  facet_grid2(rows = vars(metric), cols = vars(scenario), scales = "free_y",
              independent = "y", switch = "y") +
  labs(x = "Pooling strategy", colour = NULL, shape = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.placement = "outside",
        legend.position = "bottom", legend.box = "vertical")

save_fig("miss_cts_ci_all.png", fig_path, width = 30, height = 24, plot = fig)
