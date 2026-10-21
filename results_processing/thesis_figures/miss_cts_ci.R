##########
# title: figures for the thesis chapter - missing covariates, MI intervals
##########
# Thesis figure for the MI confidence-interval example (missing/ci_example/,
# continuous, MAR, n = 500): marginal coverage, simultaneous coverage and
# interval length by combine_mi_ci() pooling strategy, scenarios 1-4, on all,
# complete and incomplete units.
#
# - An incomplete unit's interval targets tau at its unamputed covariates, which
#   the data cannot pin down.
# - Simultaneous coverage is 0/1 per run (binomial MCSE); the others take the
#   general MCSE.
# - Every panel has its own y scale; coverage panels include the nominal 0.95.
# - Estimators keep their colours from the other figures, though only four have
#   intervals.
#
# Writes ../results/thesis_figures/miss_cts_ci/miss_cts_ci_all.png

library(here)
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
# models_figure(). The metric names are the row strips, switched to the left.
fig <- point_range_plot(
  metrics_rows,
  "value",
  NULL,
  x = "strategy",
  colour = "model",
  shape = "units",
  shape_palette = scale_shape_manual(values = UNIT_SHAPES),
  facet_rows = "metric",
  facet_scales = "free_y",
  independent_y = TRUE,
  facet_switch = "y",
  hline = nominal,
  palette = scale_colour_manual(values = pal_values, limits = model_levels),
  dodge_width = 0.85,
  point_size = 1.5,
  ci_alpha = 0.5
) +
  labs(x = "Pooling strategy", colour = NULL, shape = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.placement = "outside",
        legend.position = "bottom", legend.box = "vertical")

save_fig("miss_cts_ci_all.png", fig_path, width = 30, height = 24, plot = fig)
