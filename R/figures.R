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
  causal_forest         = "Causal forest",
  dr_random_forest      = "DR-RandomForest",
  dr_oracle             = "DR-oracle",
  dr_semi_oracle        = "DR-semi-oracle",
  dr_superlearner       = "DR-SuperLearner",
  causal_forest_inbuilt = "Causal forest (inbuilt CI)"
)

METHOD_LABELS <- c(
  complete_cases      = "Complete case",
  mean_imputation     = "Single mean imputation",
  missforest          = "Single forest-based imputation",
  regression          = "Single model-based imputation",
  missing_indicator   = "Missing indicators",
  IPW                 = "IPW",
  multiple_imputation = "Multiple imputation (rf)",
  none                = "Inbuilt missingness handling",
  complete_data       = "Complete data (reference)"
)

# The DGMs now call these MNAR / MNAR-Y throughout. AUX / AUX-Y are the older
# spelling, still present in results generated before the rename, so both map to
# the same display label.
MECHANISM_LABELS <- c(
  MAR       = "MAR",
  MNAR      = "MNAR",
  `MNAR-Y`  = "MNAR-Y",
  AUX       = "MNAR",
  `AUX-Y`   = "MNAR-Y"
)

# the four scenarios the sample-size chapters report, by their scenario index
SS_SCENARIO_LABELS <- c(`1` = "Null", `3` = "Simple", `8` = "Complex",
                        `9` = "Non-linear")

#' Recode a column to its display labels and make it an ordered factor
#'
#' Levels follow the order of `labels`, so the legend order is controlled in one
#' place. Values with no entry in `labels` are left as they are.
label_factor <- function(x, labels) {
  chr <- as.character(x)
  mapped <- ifelse(chr %in% names(labels), labels[chr], chr)
  factor(mapped, levels = unique(c(labels[names(labels) %in% chr],
                                   setdiff(mapped, labels))))
}

#' Apply the standard label recoding to a metrics tibble
#'
#' Only touches the columns that are present, so it works for both the
#' sample-size studies (scenario, n, model) and the missing-data ones (which add
#' mechanism and method).
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
  if (!is.null(scenario_labels) && "scenario" %in% names(metrics)) {
    metrics$scenario <- label_factor(metrics$scenario, scenario_labels)
  }
  metrics
}

# ---- summaries --------------------------------------------------------------

#' Mean and Monte Carlo standard error for each metric, by group
#'
#' The MCSE is sd / sqrt(number of runs), so it describes how precisely the
#' simulation has pinned down the mean - not the spread across runs.
#'
#' @param metrics one row per run
#' @param group_cols character vector of grouping columns
#' @param cols named vector: names are the output stems, values the source
#'   columns. The default keeps the stems the existing figure code already uses -
#'   note `BLP` rather than `BLP_p`, giving mean_BLP / mcse_BLP.
summarise_metrics <- function(metrics, group_cols,
                              cols = c(bias = "bias", mse = "mse", corr = "corr",
                                       BLP = "BLP_p", indep_cate = "indep_cate",
                                       indep_po = "indep_po")) {
  cols <- cols[cols %in% names(metrics)]

  out <- metrics %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(across(all_of(unname(cols)),
                     list(mean = ~ mean(.x, na.rm = TRUE),
                          mcse = ~ sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))),
                     .names = "{.fn}__{.col}"),
              .groups = "drop")

  # rename mean__BLP_p -> mean_BLP etc.
  for (stem in names(cols)) {
    for (fn in c("mean", "mcse")) {
      from <- paste0(fn, "__", cols[[stem]])
      if (from %in% names(out)) names(out)[names(out) == from] <- paste0(fn, "_", stem)
    }
  }
  out
}

# ---- presentation -----------------------------------------------------------

# the palette and theme every figure in the thesis uses
drf_scale <- function() scale_colour_paletteer_d("rcartocolor::Safe")

drf_theme <- function(rotate_x = FALSE) {
  th <- theme_minimal()
  if (rotate_x) th <- th + theme(axis.text.x = element_text(angle = 90))
  th
}

#' Save at the standard thesis figure size
save_fig <- function(filename, path, width = 21, height = 15) {
  ggsave(filename, path = path, width = width, height = height, units = "cm")
}

#' Point-and-errorbar plot of a summarised metric
#'
#' The shape both missing-data scripts use: estimator on x, one colour per
#' missing-data method, faceted by mechanism and scenario.
#'
#' @param summary output of summarise_metrics()
#' @param metric metric stem, e.g. "bias" - expects mean_<metric>/mcse_<metric>
#' @param y_lab axis label
#' @param facet_scales passed to facet_grid
point_range_plot <- function(summary, metric, y_lab, x = "model",
                             colour = "method", facet_scales = "free",
                             blank_x = FALSE) {
  mean_col <- paste0("mean_", metric)
  mcse_col <- paste0("mcse_", metric)

  p <- ggplot(summary, aes(x = .data[[x]], y = .data[[mean_col]],
                           colour = .data[[colour]],
                           ymin = .data[[mean_col]] - .data[[mcse_col]],
                           ymax = .data[[mean_col]] + .data[[mcse_col]])) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbar(position = position_dodge(width = 0.5), linewidth = 0.3) +
    facet_grid(rows = vars(mechanism), cols = vars(scenario),
               scales = facet_scales) +
    drf_scale() +
    drf_theme() +
    labs(y = y_lab, x = "Model")

  p + if (blank_x) theme(axis.text.x = element_blank()) else
    theme(axis.text.x = element_text(angle = 90))
}

#' Paired bias and MSE panels for one scenario
#'
#' Was defined identically in miss_cts.R and miss_bin.R. Free axes per scenario,
#' which is the point: the shared-scale versions compress the differences.
make_bm_plots <- function(summary, scen) {
  keep <- filter(summary, scenario == scen)
  list(
    b_plot = point_range_plot(keep, "bias", "Bias", blank_x = TRUE),
    m_plot = point_range_plot(keep, "mse", "MSE", facet_scales = "free_y",
                              blank_x = TRUE) +
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
