##########
# title: figures for the thesis chapter - model evaluation
##########
# Thesis figures for the model-evaluation strategies tree: 9 candidate CATE
# learners per run, picked by each proxy's argmin. Proxies are DR risk and the
# IF score (inflt), with estimated or fixed (0.5) propensity, x 3 fold regimes
# x 2 nuisance pipelines; calibration scores are left out. Selection and
# references come from model_evaluation/me_selection.R, shared with
# thesis_tables/me_tables.R.
#
# - Hit and regret are per (scenario, n, run, proxy), shown as the mean with a
#   95% Monte Carlo interval (binomial MCSE for hit). The best fixed choice is
#   picked with hindsight, so it is an optimistic bound on any fixed rule.
# - The rank figures' denominator is every run with a true PEHE: an unranked
#   run is a grey share rather than dropping out, so grey means a failure.
# - The axis figure averages within run first, so its MCSE counts runs.
#
# Writes to ../results/thesis_figures/model_eval/: me_top1, me_regret,
# me_axis_regret, me_rankings, me_candidate_ranks (.png).
#
#   Rscript model_eval.R [metrics.RDS] [figure dir]   # mock smoke run

library(here)
source(here("R", "figures.R"))
source(here("model_evaluation", "me_selection.R"))

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

# display labels for me_selection.R's score types, fold regimes and pipelines.
# The order of each vector is the legend order.
SCORE_LABELS <- c(dr = "DR", dr05 = "DR-π̃-fixed", inflt = "IF",
                  inflt05 = "IF-π̃-fixed")
# a light and a dark shade per score family, the dark one fixed-π. All four
# are rcartocolor::Safe colours, the thesis palette. Names are the labels above.
METRIC_COLOURS <- setNames(c("#88CCEE", "#332288", "#CC6677", "#882255"),
                           SCORE_LABELS)
FOLD_LABELS <- c(cv_shared = "CV", holdout = "Holdout", whole = "Whole")
PIPELINE_LABELS <- c(xgb = "XGBoost", automl = "AutoML")
PIPELINE_SHAPES <- setNames(c(16, 17), PIPELINE_LABELS)

# the design-axis figure's axes: the column each is read from and its levels
AXES <- list(
  Metric = list(col = "score_family", labels = c(dr = "DR", inflt = "IF")),
  Propensity = list(col = "pi_regime",
                    labels = c(estimated = "π̃ estimated",
                               fixed = "π̃ fixed")),
  `Data splitting` = list(col = "fold_type", labels = FOLD_LABELS),
  Machine = list(col = "pipeline", labels = PIPELINE_LABELS)
)

REF_LINETYPES <- c(`Random pick` = "solid", `Best fixed choice` = "dashed")

FIG_HEIGHT <- 18

# ---- data --------------------------------------------------------------------

metrics <- me_read_metrics(metrics_file, SS_SCENARIO_LABELS)

# one row per (scenario, n, run, proxy): hit, regret and the pick's true rank
selection <- me_selection(metrics) %>%
  mutate(
    metric = label_factor(score_type, SCORE_LABELS),
    fold = label_factor(fold_type, FOLD_LABELS),
    machine = label_factor(pipeline, PIPELINE_LABELS)
  )

# ---- regret references -------------------------------------------------------

refs <- me_references(metrics) %>%
  group_by(scenario, n, ref) %>%
  summarise(value = mean(regret), .groups = "drop") %>%
  mutate(ref = factor(c(random = "Random pick",
                        best_fixed = "Best fixed choice")[as.character(ref)],
                      levels = names(REF_LINETYPES))) %>%
  arrange(ref)

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

#' One proxy metric: n on x, score type by colour (METRIC_COLOURS), pipeline by
#' shape, joined across n, fold regime down, scenario across
#'
#' @param metric "hit" or "regret"
#' @param y_lab axis label
#' @param hline reference line, as point_range_plot()
proxy_figure <- function(metric, y_lab, hline = NULL) {
  point_range_plot(
    proxy_summary, metric, y_lab,
    x = "n", colour = "metric", shape = "machine",
    facet_rows = "fold", facet_cols = "scenario", facet_scales = "free_y",
    palette = scale_colour_manual(values = METRIC_COLOURS),
    shape_palette = scale_shape_manual(values = PIPELINE_SHAPES),
    hline = hline, dodge_width = 0.8, point_size = 1.5
  ) +
    labs(x = "Sample size", colour = "Metric", shape = "Machine") +
    theme(legend.position = "bottom", legend.box = "vertical",
          axis.text.x = element_text(angle = 0))
}

top1_fig <- proxy_figure("hit", "Top-1 selection accuracy",
                         hline = 1 / ME_N_CANDIDATES)
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

axis_runs <- me_axis_runs(selection, AXES, "regret") %>%
  mutate(n = n_label(n))

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

# ---- true rank of the pick ---------------------------------------------------

# fixed-π proxies only, so the propensity is left out of the bar labels. Bars
# are grouped by score, then fold regime, then pipeline.
rank_levels <- with(
  expand.grid(p = PIPELINE_LABELS, f = FOLD_LABELS, s = c("DR", "IF"),
              stringsAsFactors = FALSE),
  paste(s, f, p, sep = " / ")
)

NOT_RANKED <- "Not ranked"

# the runs with a true PEHE, per (scenario, n): each bar's denominator
total_runs <- me_candidates(metrics) %>%
  distinct(scenario, n, run) %>%
  count(scenario, n, name = "total")

ranked <- selection %>%
  filter(pi_regime == "fixed") %>%
  mutate(
    proxy = factor(
      paste(c(dr = "DR", inflt = "IF")[score_family], FOLD_LABELS[fold_type],
            PIPELINE_LABELS[pipeline], sep = " / "),
      levels = rank_levels
    ),
    rank = as.character(pick_rank)
  ) %>%
  count(scenario, n, proxy, rank, name = "runs")

# over every proxy, so one that failed in every run is a whole grey bar
not_ranked <- tidyr::expand_grid(
  total_runs, proxy = factor(rank_levels, levels = rank_levels)
) %>%
  left_join(count(ranked, scenario, n, proxy, wt = runs, name = "ranked"),
            by = c("scenario", "n", "proxy")) %>%
  transmute(scenario, n, proxy, rank = NOT_RANKED,
            runs = total - coalesce(ranked, 0L))

# every possible rank, so a rank's colour does not depend on which ranks occur
rank_names <- as.character(seq_len(ME_N_CANDIDATES))
rank_dist <- bind_rows(ranked, not_ranked) %>%
  filter(runs > 0) %>%
  left_join(total_runs, by = c("scenario", "n")) %>%
  mutate(
    prop = runs / total,
    rank = factor(rank, levels = c(rank_names, NOT_RANKED)),
    n = n_label(n)
  )

unranked <- filter(rank_dist, rank == NOT_RANKED)
if (nrow(unranked) > 0) {
  message("bars with unranked runs (", nrow(unranked), "):")
  print(as.data.frame(select(unranked, scenario, n, proxy, runs, total)),
        row.names = FALSE)
}

# viridis rather than the thesis palette: the rank is ordinal, and Safe would
# give it unordered hues. Reversed, so the best ranks are the bright ones. Not
# ranked is the last level, so it stacks at the foot of the bar, below rank 9.
rank_colours <- c(
  setNames(scales::viridis_pal(direction = -1)(length(rank_names)),
           rank_names),
  setNames("grey75", NOT_RANKED)
)

# show.legend = TRUE, or a rank no pick took gets an empty legend key
rank_fig <- ggplot(rank_dist, aes(x = proxy, y = prop, fill = rank)) +
  geom_col(show.legend = TRUE) +
  facet_grid(rows = vars(n), cols = vars(scenario)) +
  scale_fill_manual(values = rank_colours, drop = FALSE) +
  drf_theme() +
  labs(x = "Proxy score", y = "Proportion of runs",
       fill = "True rank of pick") +
  guides(fill = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
save_fig("me_rankings.png", fig_path, height = FIG_HEIGHT, plot = rank_fig)

# ---- true rank of each candidate ---------------------------------------------

# me_config.R's CANDIDATE_MODELS, in display order. Not sourced from there:
# me_config.R sets up the study run.
CANDIDATE_LEVELS <- c("rf1", "rf2", "rf3", "net1", "net2", "net3",
                      "SL1", "SL2", "SL3")

cand_ranked <- me_candidates(metrics) %>%
  mutate(model = factor(model, levels = CANDIDATE_LEVELS),
         rank = as.character(pick_rank)) %>%
  count(scenario, n, model, rank, name = "runs")

# over every candidate, as not_ranked above, with the same total_runs
cand_not_ranked <- tidyr::expand_grid(
  total_runs, model = factor(CANDIDATE_LEVELS, levels = CANDIDATE_LEVELS)
) %>%
  left_join(count(cand_ranked, scenario, n, model, wt = runs, name = "ranked"),
            by = c("scenario", "n", "model")) %>%
  transmute(scenario, n, model, rank = NOT_RANKED,
            runs = total - coalesce(ranked, 0L))

cand_dist <- bind_rows(cand_ranked, cand_not_ranked) %>%
  filter(runs > 0) %>%
  left_join(total_runs, by = c("scenario", "n")) %>%
  mutate(
    prop = runs / total,
    rank = factor(rank, levels = c(rank_names, NOT_RANKED)),
    n = n_label(n)
  )

cand_fig <- ggplot(cand_dist, aes(x = model, y = prop, fill = rank)) +
  geom_col(show.legend = TRUE) +
  facet_grid(rows = vars(n), cols = vars(scenario)) +
  scale_fill_manual(values = rank_colours, drop = FALSE) +
  drf_theme() +
  labs(x = "Candidate", y = "Proportion of runs", fill = "True rank") +
  guides(fill = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
save_fig("me_candidate_ranks.png", fig_path, height = FIG_HEIGHT,
         plot = cand_fig)
