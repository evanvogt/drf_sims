##########
# title: shared plotting helpers for the thesis figures
##########
# results_processing/thesis_figures/ repeated itself: cts_ss.R and bin_ss.R
# differed in 18 of 280 lines once the outcome prefix was normalised, miss_cts.R
# and miss_bin.R in 64 of 140, and make_bm_plots() was defined identically in
# both of the latter.
#
# What is shared is the presentation, not the analysis: the display labels, the
# palette and theme, the mean +/- Monte Carlo standard error summary, and the
# figure size. Each study script keeps its own paths, filters and plot choices.

# only what the helpers below actually need. Scripts that use patchwork (`/`,
# plot_layout), ggridges, scales or purrr load those themselves.
library(dplyr)
library(tidyr)
library(ggplot2)
library(paletteer)

# ---- display labels ---------------------------------------------------------
# One place to rename an estimator for the write-up.

MODEL_LABELS <- c(
  causal_forest = "Causal forest",
  dr_random_forest = "DR-RandomForest",
  dr_oracle = "DR-oracle",
  dr_semi_oracle = "DR-semi-oracle",
  dr_superlearner = "DR-SuperLearner",
  causal_forest_inbuilt = "Causal forest (inbuilt CI)"
)

METHOD_LABELS <- c(
  complete_cases = "Complete case",
  mean_imputation = "Single mean imputation",
  missforest = "Single forest-based imputation",
  regression = "Single model-based imputation",
  missing_indicator = "Missing indicators",
  IPW = "IPW",
  multiple_imputation = "Multiple imputation (rf)",
  none = "Inbuilt missingness handling",
  complete_data = "Complete data (reference)"
)

# The DGMs now call these MNAR / MNAR-Y throughout. AUX / AUX-Y are the older
# spelling, still present in results generated before the rename, so both map to
# the same display label.
MECHANISM_LABELS <- c(
  MAR = "MAR",
  MNAR = "MNAR",
  `MNAR-Y` = "MNAR-Y",
  AUX = "MNAR",
  `AUX-Y` = "MNAR-Y"
)

# the four scenarios the sample-size chapters report, by their scenario index
SS_SCENARIO_LABELS <- c(
  `1` = "Null",
  `3` = "Simple",
  `8` = "Complex",
  `9` = "Non-linear"
)

# the missing-data studies' renumbered five-scenario set. Was copy-pasted
# verbatim into all three thesis_figures/miss_*.R scripts, minus scenario 3.
# continuous/ and binary/ never run 3 (their grid is 1, 2, 4, 5), so it drops out
# of those by droplevels(); ci_example runs all five and keeps it.
MISS_SCENARIO_LABELS <- c(
  `1` = "Null",
  `2` = "Simple",
  `3` = "Two HTE vars",
  `4` = "Complex",
  `5` = "Non-linear"
)

# the three ways combine_mi_ci() (R/bootstrap_ci.R) pools per-imputation
# bootstrap intervals. These are the values cts_miss_ci_metrics.R writes to the
# `strategy` column.
STRATEGY_LABELS <- c(
  pooled = "Pooled quantiles",
  mib = "MI bootstrap (Rubin)",
  hybrid = "Hybrid"
)

#' Recode a column to its display labels and make it an ordered factor
#'
#' Levels follow the order of `labels`, so the legend order is controlled in one
#' place. Values with no entry in `labels` are left as they are.
label_factor <- function(x, labels) {
  chr <- as.character(x)
  mapped <- ifelse(chr %in% names(labels), labels[chr], chr)
  factor(
    mapped,
    levels = unique(c(labels[names(labels) %in% chr], setdiff(mapped, labels)))
  )
}

#' Apply the standard label recoding to a metrics tibble
#'
#' Only touches the columns that are present, so it works for the sample-size
#' studies (scenario, n, model), the missing-data ones (which add mechanism and
#' method) and the MI confidence-interval one (which adds strategy).
apply_labels <- function(metrics, scenario_labels = NULL) {
  if ("model" %in% names(metrics)) {
    metrics$model <- label_factor(metrics$model, MODEL_LABELS)
  }
  if ("method" %in% names(metrics)) {
    metrics$method <- label_factor(metrics$method, METHOD_LABELS)
  }
  if ("mechanism" %in% names(metrics)) {
    metrics$mechanism <- label_factor(metrics$mechanism, MECHANISM_LABELS)
  }
  if ("strategy" %in% names(metrics)) {
    metrics$strategy <- label_factor(metrics$strategy, STRATEGY_LABELS)
  }
  if (!is.null(scenario_labels) && "scenario" %in% names(metrics)) {
    metrics$scenario <- label_factor(metrics$scenario, scenario_labels)
  }
  metrics
}

# ---- summaries --------------------------------------------------------------

#' Mean and Monte Carlo standard error for each metric, by group
#'
#' For a continuous measure (bias, MSE, a p-value, ...) the MCSE is
#' sd(x) / sqrt(number of non-NA runs) - it describes how precisely the
#' simulation has pinned down the mean, not the spread across runs. For a
#' genuine 0/1 indicator measured once per run (coverage of the whole band,
#' a hit/miss pick, ...) the correct MCSE is instead the binomial
#' sqrt(p(1-p)/n) - using the general sd-based formula there is not wrong in
#' the limit, but doesn't match the convention `rsimsum`/Morris, White &
#' Crowther (2019) document for these measures. List any such column's stem
#' in `binomial` to get that treatment instead. Do NOT list a column that is
#' a proportion averaged over many units *within* a run (e.g. marginal
#' coverage, sign accuracy) - those are not Bernoulli at the run level, and
#' the binomial formula would badly overstate their MCSE.
#'
#' @param metrics one row per run
#' @param group_cols character vector of grouping columns
#' @param cols named vector: names are the output stems, values the source
#'   columns. The default keeps the stems the existing figure code already uses -
#'   note `BLP` rather than `BLP_p`, giving mean_BLP / mcse_BLP.
#' @param binomial character vector of stems (from `names(cols)`) that are
#'   genuine per-run 0/1 indicators and should get the binomial MCSE formula
#'   instead of the general one. Empty by default.
summarise_metrics <- function(
  metrics,
  group_cols,
  cols = c(
    bias = "bias",
    mse = "mse",
    corr = "corr",
    BLP = "BLP_p",
    indep_cate = "indep_cate",
    indep_po = "indep_po",
    rel_ate_bias = "rel_ate_bias",
    rel_bias_cate = "rel_bias_cate"
  ),
  binomial = character()
) {
  cols <- cols[cols %in% names(metrics)]
  cont_cols <- cols[!names(cols) %in% binomial]
  binom_cols <- cols[names(cols) %in% binomial]

  out <- metrics %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      across(
        all_of(unname(cont_cols)),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          mcse = ~ sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
        ),
        .names = "{.fn}__{.col}"
      ),
      across(
        all_of(unname(binom_cols)),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          mcse = ~ {
            p <- mean(.x, na.rm = TRUE)
            sqrt(p * (1 - p) / sum(!is.na(.x)))
          }
        ),
        .names = "{.fn}__{.col}"
      ),
      .groups = "drop"
    )

  # rename mean__BLP_p -> mean_BLP etc.
  for (stem in names(cols)) {
    for (fn in c("mean", "mcse")) {
      from <- paste0(fn, "__", cols[[stem]])
      if (from %in% names(out)) {
        names(out)[names(out) == from] <- paste0(fn, "_", stem)
      }
    }
  }
  out
}

# ---- presentation -----------------------------------------------------------

# the palette and theme every figure in the thesis uses
drf_scale <- function() scale_colour_paletteer_d("rcartocolor::Safe")

drf_theme <- function(rotate_x = FALSE) {
  th <- theme_light() +
    theme(
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(colour = "black")
    )
  if (rotate_x) {
    th <- th + theme(axis.text.x = element_text(angle = 90))
  }
  th
}

#' Save at the standard thesis figure size
#'
#' @param plot defaults to ggplot2::last_plot(), same as a bare ggsave() call
#'   - fine when the plot was just built immediately before this call, but
#'   pass it explicitly when saving several plots pulled out of a list (e.g.
#'   poster_figs/explore_gallery.R), where "last plot" is fragile
save_fig <- function(
  filename,
  path,
  width = 21,
  height = 15,
  plot = ggplot2::last_plot()
) {
  ggsave(
    filename,
    plot = plot,
    path = path,
    width = width,
    height = height,
    units = "cm"
  )
}

#' Point-and-errorbar plot of a summarised metric
#'
#' The shape both missing-data scripts use: estimator on x, one colour per
#' missing-data method, faceted by mechanism and scenario.
#'
#' @param summary output of summarise_metrics()
#' @param metric metric stem, e.g. "bias" - expects mean_<metric>/mcse_<metric>
#' @param y_lab axis label
#' @param facet_rows,facet_cols character vectors of column names to facet by
#'   (NULL to skip that dimension), passed to facet_grid() via vars() - as
#'   lollipop_plot(). Default to the mechanism/scenario grid this function
#'   originally hardcoded, so existing callers are unaffected.
#' @param facet_scales passed to facet_grid
#' @param alpha error bar is a (1 - alpha) CI, i.e. mean +/- qnorm(1 -
#'   alpha/2) x MCSE - not a raw +/- 1x MCSE, which is only a ~68% interval
#'   (see continuous/cts_results.R's summary_plot() for the verification of
#'   why the raw form is wrong). Not to be confused with `ci_alpha` below -
#'   this `alpha` is a significance level, not a plotting transparency.
#' @param line add a geom_line connecting each colour group's points across x
#'   (drawn under the points) - only meaningful when x has a natural order
#'   (e.g. sample size), not a bare category. When `shape` is also set, each
#'   (colour, shape) pair gets its own line rather than one line per colour
#'   that jumps between shapes.
#' @param ci_alpha opacity of the error bar only (points and the line, if
#'   drawn, stay fully opaque) - useful to de-emphasise the CI in a busy,
#'   many-colour panel
#' @param palette the colour scale layer; defaults to the thesis-standard
#'   drf_scale(), override with e.g. scale_colour_manual(values = ...) for a
#'   study-specific palette
#' @param shape optional column name mapped to point/line shape, alongside
#'   `colour` - e.g. distinguishing two evaluation kinds within one colour
#'   group. NULL (the default) omits the aesthetic entirely, so every
#'   existing caller is unaffected.
#' @param shape_palette optional scale layer for `shape` (e.g.
#'   scale_shape_manual(values = ...)); ignored when `shape` is NULL
#' @param hline reference line, dashed - 0 for bias-type metrics (the
#'   original hardcoded value), but e.g. 0.95 for a nominal-coverage panel
point_range_plot <- function(
  summary,
  metric,
  y_lab,
  x = "model",
  colour = "method",
  facet_rows = "mechanism",
  facet_cols = "scenario",
  facet_scales = "free",
  blank_x = FALSE,
  alpha = 0.05,
  line = FALSE,
  ci_alpha = 1,
  palette = drf_scale(),
  shape = NULL,
  shape_palette = NULL,
  hline = 0
) {
  mean_col <- paste0("mean_", metric)
  mcse_col <- paste0("mcse_", metric)
  z <- qnorm(1 - alpha / 2)
  x_lab <- if (x == "model") "Model" else x

  p <- ggplot(
    summary,
    aes(
      x = .data[[x]],
      y = .data[[mean_col]],
      colour = .data[[colour]],
      ymin = .data[[mean_col]] - z * .data[[mcse_col]],
      ymax = .data[[mean_col]] + z * .data[[mcse_col]]
    )
  )

  if (!is.null(shape)) {
    p <- p + aes(shape = .data[[shape]])
  }

  p <- p + geom_hline(yintercept = hline, linetype = "dashed")

  if (line) {
    p <- p +
      geom_line(
        aes(
          group = if (is.null(shape)) {
            .data[[colour]]
          } else {
            interaction(.data[[colour]], .data[[shape]])
          }
        ),
        position = position_dodge(width = 0.5),
        linewidth = 0.5,
        alpha = 0.5
      )
  }

  p <- p +
    geom_point(position = position_dodge(width = 0.5), size = 1.5) +
    geom_errorbar(
      position = position_dodge(width = 0.5),
      linewidth = 0.3,
      alpha = ci_alpha
    ) +
    facet_grid(
      rows = vars(!!!rlang::syms(facet_rows)),
      cols = vars(!!!rlang::syms(facet_cols)),
      scales = facet_scales
    ) +
    palette +
    drf_theme() +
    labs(y = y_lab, x = x_lab)

  if (!is.null(shape_palette)) {
    p <- p + shape_palette
  }

  p +
    if (blank_x) {
      theme(axis.text.x = element_blank())
    } else {
      theme(axis.text.x = element_text(angle = 90))
    }
}

#' Lollipop plot of a summarised metric: stick + point + a narrow MCSE bracket
#'
#' The rsimsum-style alternative to point_range_plot(): a segment from a
#' reference value (`hline`) up to the point estimate, a point at the tip, and
#' a narrow errorbar bracket at the tip for the (1 - alpha) MCSE-based CI -
#' rather than one wide errorbar spanning the full range. Useful whenever the
#' story is deviation from a target (0 bias, 1 correlation, a target
#' power/coverage), which several of the ICTMC poster panels want - see
#' results_processing/poster_figs/R/candidates.R.
#'
#' Facet rows/cols are parameters rather than hardcoded, unlike
#' point_range_plot()'s fixed mechanism/scenario grid, so the same helper
#' covers both that shape and others (e.g. the sample-size studies' facet by
#' scenario only, with sample size on the x axis).
#'
#' @param summary output of summarise_metrics()
#' @param metric metric stem, e.g. "bias" - expects mean_<metric>/mcse_<metric>
#' @param y_lab axis label
#' @param hline reference value the stick starts from (0 for bias-type
#'   metrics, 1 for correlation-type metrics where 1 is ideal, a target
#'   power/coverage for test/CI panels)
#' @param x,colour aesthetics, as in point_range_plot()
#' @param facet_rows,facet_cols character vectors of column names to facet by
#'   (NULL to skip that dimension), passed to facet_grid() via vars()
#' @param facet_scales passed to facet_grid
#' @param alpha the bracket is a (1 - alpha) CI, i.e. mean +/- qnorm(1 -
#'   alpha/2) x MCSE - see point_range_plot() for why not a raw +/- 1x MCSE
#' @param bracket_width width of the MCSE bracket, narrow relative to the
#'   dodge width so it reads as a cap on the stick, not a second errorbar
lollipop_plot <- function(
  summary,
  metric,
  y_lab,
  hline = 0,
  x = "model",
  colour = "method",
  facet_rows = "mechanism",
  facet_cols = "scenario",
  facet_scales = "free",
  blank_x = FALSE,
  alpha = 0.05,
  bracket_width = 0.15
) {
  mean_col <- paste0("mean_", metric)
  mcse_col <- paste0("mcse_", metric)
  z <- qnorm(1 - alpha / 2)

  summary <- summary %>%
    mutate(
      .mean = .data[[mean_col]],
      .lo = .data[[mean_col]] - z * .data[[mcse_col]],
      .hi = .data[[mean_col]] + z * .data[[mcse_col]]
    )

  p <- ggplot(summary, aes(x = .data[[x]], colour = .data[[colour]])) +
    geom_hline(yintercept = hline, linetype = "dashed") +
    # the stick: geom_segment() + position_dodge() dodges x but leaves xend
    # untouched, so every group's stick converges on one undodged point -
    # geom_linerange() has no xend, just a vertical range at one (dodged) x,
    # so it dodges the same clean way the errorbar bracket below it does
    geom_linerange(
      aes(ymin = pmin(hline, .mean), ymax = pmax(hline, .mean)),
      position = position_dodge(width = 0.5),
      linewidth = 0.5
    ) +
    geom_errorbar(
      aes(ymin = .lo, ymax = .hi),
      position = position_dodge(width = 0.5),
      width = bracket_width,
      linewidth = 0.3,
      alpha = 0.5
    ) +
    geom_point(
      aes(y = .mean),
      position = position_dodge(width = 0.5),
      size = 2
    ) +
    facet_grid(
      rows = vars(!!!rlang::syms(facet_rows)),
      cols = vars(!!!rlang::syms(facet_cols)),
      scales = facet_scales
    ) +
    drf_scale() +
    drf_theme() +
    labs(y = y_lab, x = x)

  p +
    if (blank_x) {
      theme(axis.text.x = element_blank())
    } else {
      theme(axis.text.x = element_text(angle = 90))
    }
}

#' Paired bias and MSE panels for one scenario
#'
#' Was defined identically in miss_cts.R and miss_bin.R. Free axes per scenario,
#' which is the point: the shared-scale versions compress the differences.
make_bm_plots <- function(summary, scen) {
  keep <- filter(summary, scenario == scen)
  list(
    b_plot = point_range_plot(keep, "bias", "Bias", blank_x = TRUE),
    m_plot = point_range_plot(
      keep,
      "mse",
      "MSE",
      facet_scales = "free_y",
      blank_x = TRUE
    ) +
      theme(axis.text = element_blank())
  )
}

#' Boxplot of a per-run metric against sample size
#'
#' The shape both sample-size scripts use.
distribution_plot <- function(metrics, metric, y_lab, title = NULL) {
  ggplot(metrics, aes(x = n, y = .data[[metric]], colour = model)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_boxplot(fill = "transparent", outlier.shape = NA) +
    facet_wrap(~scenario) +
    drf_scale() +
    drf_theme() +
    labs(title = title, y = y_lab, x = "Sample size")
}
