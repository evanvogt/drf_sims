##########
# title: figures for the thesis chapter - model evaluation
##########
# The model-evaluation study (model_evaluation/), the strategies tree: 9
# candidate CATE learners per run, each scored by proxy losses built from an
# independent nuisance pipeline, and selected by the proxy's argmin. The
# exploratory version of these figures is model_evaluation/me_results.qmd; the
# selection and baseline derivations here are the same.
#
# Proxies: the four PEHE proxies, DR risk and the influence-corrected score
# (IF), each with the pipeline's estimated propensity or with it fixed at the
# true 0.5 (the `05` score types, model_evaluation/me_metrics.R), crossed with
# the three fold regimes (CV = cv_shared, holdout, whole) and the two nuisance
# pipelines (XGBoost, AutoML). The calibration scores are left out.
#
# Per (scenario, n, run, proxy):
#   - hit, whether the proxy's pick is the candidate with the lowest true PEHE;
#   - regret, the pick's true PEHE minus the lowest true PEHE in the run.
# Mean over runs with a 95% Monte Carlo interval (binomial MCSE for hit).
#
# Regret reference lines, per (scenario, n):
#   - random pick, the expected regret of a uniformly random candidate;
#   - best fixed choice, the single candidate with the lowest mean regret,
#     chosen with hindsight from the truth, so an optimistic bound on any
#     fixed rule.
# Top-1 accuracy's reference is random-pick accuracy, 1 / 9.
#
# The design-axis figure marginalises regret over the other axes: each level
# (DR vs IF, estimated vs fixed propensity, fold regime, pipeline) is averaged
# within run first, then over runs, so its MCSE counts runs, not runs x the
# proxies sharing that level (me_results.qmd pools the latter, which
# understates it).
#
# Labels, palette, summaries and figure sizing come from R/figures.R. This
# script carries only the paths and this study's filters.
#
# Writes to ../results/thesis_figures/model_eval/:
#   me_top1.png, me_regret.png, me_axis_regret.png
#
# Optional arguments, for a smoke run on mock metrics:
#   Rscript model_eval.R [metrics.RDS] [figure dir]

library(here)
source(here("R", "figures.R"))

# paths
args <- commandArgs(trailingOnly = TRUE)
path <- here()
metrics_file <- if (length(args) >= 1) {
  args[1]
} else {
  file.path(dirname(path), "results", "correlated",
            "model_evaluation_strategies", "me_strat_metrics.RDS")
}
fig_path <- if (length(args) >= 2) {
  args[2]
} else {
  file.path(dirname(path), "results", "thesis_figures", "model_eval")
}
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

KEY_COLS <- c("scenario", "n", "run", "model")
N_CANDIDATES <- 9

# the score columns are <score_type>_<fold_type>_<pipeline>
# (model_evaluation/me_metrics.R). The order of each vector is the legend order.
SCORE_LABELS <- c(dr = "DR", dr05 = "DR-π-fixed", infl = "IF",
                  infl05 = "IF-π-fixed")
FOLD_LABELS <- c(cv_shared = "CV", holdout = "Holdout", whole = "Whole")
PIPELINE_LABELS <- c(xgb = "XGBoost", automl = "AutoML")
PIPELINE_SHAPES <- setNames(c(16, 17), PIPELINE_LABELS)

# the design-axis figure's axes: the column each is read from and its levels
AXES <- list(
  Score = list(col = "score_family", labels = c(dr = "DR", infl = "IF")),
  Propensity = list(col = "pi_regime",
                    labels = c(estimated = "π estimated", fixed = "π fixed")),
  `Fold regime` = list(col = "fold_type", labels = FOLD_LABELS),
  Machine = list(col = "pipeline", labels = PIPELINE_LABELS)
)

REF_LINETYPES <- c(`Random pick` = "solid", `Best fixed choice` = "dashed")

FIG_HEIGHT <- 18

# anchored on the known tokens, so cv_shared's own underscore can't misparse
score_pattern <- paste0(
  "^(", paste(names(SCORE_LABELS), collapse = "|"), ")_",
  "(", paste(names(FOLD_LABELS), collapse = "|"), ")_",
  "(", paste(names(PIPELINE_LABELS), collapse = "|"), ")$"
)

# ---- data --------------------------------------------------------------------

metrics <- readRDS(metrics_file)

expected_cols <- with(
  expand.grid(s = names(SCORE_LABELS), f = names(FOLD_LABELS),
              p = names(PIPELINE_LABELS), stringsAsFactors = FALSE),
  paste(s, f, p, sep = "_")
)
absent <- setdiff(expected_cols, names(metrics))
if (length(absent) > 0) {
  stop(metrics_file, " has no ", paste(absent, collapse = ", "),
       ". Is it the strategies tree (Rscript me_metrics.R strategies)?",
       call. = FALSE)
}

metrics <- metrics %>%
  filter(scenario %in% names(SS_SCENARIO_LABELS)) %>%
  select(all_of(KEY_COLS), true_pehe, all_of(expected_cols)) %>%
  mutate(
    scenario = label_factor(scenario, SS_SCENARIO_LABELS),
    n = factor(n, levels = sort(unique(n)))
  )

# one row per (scenario, n, run, proxy). A candidate only enters a ranking if
# both its score and its truth are there; regret is against the best candidate
# present in that run, as in me_results.qmd.
selection <- metrics %>%
  pivot_longer(all_of(expected_cols), names_to = "proxy_id",
               values_to = "score") %>%
  filter(!is.na(score), !is.na(true_pehe)) %>%
  mutate(
    score_type = sub(score_pattern, "\\1", proxy_id),
    fold_type = sub(score_pattern, "\\2", proxy_id),
    pipeline = sub(score_pattern, "\\3", proxy_id)
  ) %>%
  group_by(scenario, n, run, score_type, fold_type, pipeline) %>%
  filter(n() >= 2) %>%
  summarise(
    hit = as.numeric(model[which.min(score)] == model[which.min(true_pehe)]),
    regret = true_pehe[which.min(score)] - min(true_pehe),
    .groups = "drop"
  ) %>%
  mutate(
    score_family = sub("05$", "", score_type),
    pi_regime = ifelse(grepl("05$", score_type), "fixed", "estimated"),
    metric = label_factor(score_type, SCORE_LABELS),
    fold = label_factor(fold_type, FOLD_LABELS),
    machine = label_factor(pipeline, PIPELINE_LABELS)
  )

# ---- regret references -------------------------------------------------------

# each candidate's regret in its run, the basis of both references
candidate_regret <- metrics %>%
  filter(!is.na(true_pehe)) %>%
  group_by(scenario, n, run) %>%
  mutate(regret = true_pehe - min(true_pehe)) %>%
  ungroup()

# a uniformly random pick's expected regret is the run's mean candidate regret
random_ref <- candidate_regret %>%
  group_by(scenario, n, run) %>%
  summarise(regret = mean(regret), .groups = "drop") %>%
  group_by(scenario, n) %>%
  summarise(value = mean(regret), .groups = "drop") %>%
  mutate(ref = "Random pick")

best_fixed_ref <- candidate_regret %>%
  group_by(scenario, n, model) %>%
  summarise(value = mean(regret), .groups = "drop_last") %>%
  slice_min(value, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(ref = "Best fixed choice")

refs <- bind_rows(random_ref, select(best_fixed_ref, -model)) %>%
  mutate(ref = factor(ref, levels = names(REF_LINETYPES)))

ref_linetype <- scale_linetype_manual(name = "Reference",
                                      values = REF_LINETYPES)

# ---- top-1 accuracy and regret, per proxy -------------------------------------

proxy_summary <- summarise_metrics(
  selection,
  c("scenario", "n", "metric", "fold", "machine"),
  cols = c(hit = "hit", regret = "regret"),
  binomial = "hit",
  count_na = character()
)

#' One proxy metric: n on x, score type by colour, pipeline by shape, fold
#' regime down, scenario across
#'
#' @param metric "hit" or "regret"
#' @param y_lab axis label
#' @param hline reference line, as point_range_plot()
proxy_figure <- function(metric, y_lab, hline = NULL) {
  point_range_plot(
    proxy_summary, metric, y_lab,
    x = "n", colour = "metric", shape = "machine",
    facet_rows = "fold", facet_cols = "scenario", facet_scales = "free_y",
    shape_palette = scale_shape_manual(values = PIPELINE_SHAPES),
    hline = hline, dodge_width = 0.8, point_size = 1.5
  ) +
    labs(x = "Sample size", colour = "Metric", shape = "Machine") +
    theme(legend.position = "bottom", legend.box = "vertical",
          axis.text.x = element_text(angle = 0))
}

top1_fig <- proxy_figure("hit", "Top-1 selection accuracy",
                         hline = 1 / N_CANDIDATES)
save_fig("me_top1.png", fig_path, height = FIG_HEIGHT, plot = top1_fig)

# the references vary with n, which is on x, so each is a segment across its n
# slot rather than a full-width line; refs has no fold column, so the same
# segments repeat in every row
regret_fig <- proxy_figure("regret", "Mean regret") +
  geom_errorbar(
    data = refs,
    aes(x = n, ymin = value, ymax = value, linetype = ref),
    width = 0.9, colour = "grey30", linewidth = 0.4,
    inherit.aes = FALSE, key_glyph = "path"
  ) +
  ref_linetype
save_fig("me_regret.png", fig_path, height = FIG_HEIGHT, plot = regret_fig)

# ---- regret by design axis ---------------------------------------------------

n_label <- function(n) factor(paste("n =", n), levels = paste("n =", levels(n)))

axis_runs <- bind_rows(Map(function(axis, spec) {
  selection %>%
    mutate(axis = axis, level = unname(spec$labels[.data[[spec$col]]])) %>%
    group_by(scenario, n, run, axis, level) %>%
    summarise(regret = mean(regret), .groups = "drop")
}, names(AXES), AXES)) %>%
  mutate(
    axis = factor(axis, levels = names(AXES)),
    level = factor(level,
                   levels = unlist(lapply(AXES, function(a) unname(a$labels)))),
    n = n_label(n)
  )

axis_summary <- summarise_metrics(
  axis_runs,
  c("scenario", "n", "axis", "level"),
  cols = c(regret = "regret"),
  count_na = character()
)

axis_fig <- point_range_plot(
  axis_summary, "regret", "Mean regret",
  x = "level", colour = "axis",
  facet_rows = "n", facet_cols = "scenario", facet_scales = "free_y",
  hline = NULL, point_size = 1.5
) +
  geom_hline(
    data = mutate(refs, n = n_label(n)),
    aes(yintercept = value, linetype = ref),
    colour = "grey30", linewidth = 0.4
  ) +
  ref_linetype +
  labs(x = NULL, colour = "Design axis") +
  theme(legend.position = "bottom", legend.box = "vertical",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
save_fig("me_axis_regret.png", fig_path, height = FIG_HEIGHT, plot = axis_fig)
