##########
# title: figures for the thesis chapter - competing risks
##########
# Thesis figures for the competing-risks study: CATE bias, RMSE and Pearson,
# each at rho = 0.5, rho = 0 and their paired difference; censoring on x, model
# by colour, event by row, scenario by column. One arm per estimator family, as
# in thesis_tables/surv_tables.R's main table.
#
# - Arms target different truths (ADEMP.md, "Estimands"); shape marks the
#   net-RMST arms (ipw, csf_cs) against the RMTL-scale ones.
# - Combined (RMSTc) is in every figure but the rho = 0.5 bias and RMSE ones
#   (EVENTS_ONLY); csf_sh has none.
# - Pearson is NA where the arm's truth is constant, so those points are absent.
# - Differences are per run (shared seed across rho), MCSE sd(diff) / sqrt(pairs).
#
# Writes to ../results/thesis_figures/surv/:
#   surv_{rho05,rho0,rhodiff}_{bias,rmse,corr}.png

library(here)
source(here("R", "figures.R"))

# paths
path <- here()
res_path <- file.path(dirname(path), "results", "competing_risk")
fig_path <- file.path(dirname(path), "results", "thesis_figures", "surv")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

# the levels views, by their file-name tag; "rhodiff" is the third view
RHOS <- c(rho05 = 0.5, rho0 = 0)

TARGETS <- c("Event 1", "Event 2", "Combined")

# figures shown without the Combined row, by "<view>_<metric>"
EVENTS_ONLY <- c("rho05_bias", "rho05_rmse")

ESTIMAND_LABELS <- c(net = "Net RMST", sub = "RMTL")
ESTIMAND_SHAPES <- setNames(c(17, 16), ESTIMAND_LABELS)

# second line of a difference figure's y labels
DIFF_LAB <- "ρ = 0.5 − ρ = 0"

# seven scenario columns, so wider than save_fig()'s standard 21cm
FIG_WIDTH <- 33
FIG_HEIGHT <- 20
EVENTS_ONLY_HEIGHT <- 15  # two event rows rather than three

# per metric: the levels' and the differences' y labels, and the levels'
# reference line (the differences' is 0 throughout)
FIGURES <- list(
  bias = list(y_lab = "Bias of the CATE (days)",
              diff_lab = "Δ bias of the CATE (days)", hline = 0),
  rmse = list(y_lab = "RMSE of the CATE (days)",
              diff_lab = "Δ RMSE of the CATE (days)", hline = NULL),
  corr = list(y_lab = "Pearson correlation with the true CATE",
              diff_lab = "Δ Pearson correlation", hline = NULL)
)

# ---- data --------------------------------------------------------------------

metrics <- readRDS(file.path(res_path, "surv_metrics.RDS"))
if (!"rho" %in% names(metrics)) {
  stop("surv_metrics.RDS has no rho column: it predates the 2026-10-01 grid. ",
       "Rerun study.R collect and surv_metrics.R on the current results.",
       call. = FALSE)
}
s1_corr <- metrics$corr[metrics$scenario == 1]
if (length(s1_corr) > 0 && all(s1_corr %in% 0)) {
  stop("scenario 1's Pearson correlations are all 0, cate_metrics()'s ",
       "null-scenario placeholder: surv_metrics.RDS predates the 2026-10-06 ",
       "fix. Rerun surv_metrics.R.", call. = FALSE)
}

metrics <- metrics %>%
  filter(framework %in% names(SURV_ARM_LABELS), target %in% TARGETS)

#' One view's per-run rows: one rho's, or the paired rho differences
#'
#' @param view "rho05", "rho0" or "rhodiff"
view_runs <- function(view) {
  if (view == "rhodiff") {
    paired_rho_diff(
      metrics,
      names(FIGURES),
      keys = c("scenario", "n", "censoring", "framework", "target", "run")
    )
  } else {
    filter(metrics, rho == RHOS[[view]])
  }
}

#' Mean and MCSE of every metric by scenario, censoring, model and event
summarise_view <- function(runs) {
  runs %>%
    mutate(
      model = factor(SURV_ARM_LABELS[framework],
                     levels = unname(SURV_ARM_LABELS)),
      estimand = factor(
        if_else(framework %in% SURV_NET_ARMS,
                ESTIMAND_LABELS[["net"]], ESTIMAND_LABELS[["sub"]]),
        levels = ESTIMAND_LABELS
      ),
      scenario = label_factor(scenario, SURV_SCENARIO_PLOT_LABELS),
      censoring = label_factor(censoring, SURV_CENSORING_LABELS),
      target = factor(target, levels = TARGETS)
    ) %>%
    summarise_metrics(
      c("scenario", "censoring", "model", "estimand", "target"),
      cols = setNames(names(FIGURES), names(FIGURES)),
      count_na = character()
    )
}

# ---- figures -----------------------------------------------------------------

#' One metric for one view
#'
#' @param summary summarise_view() output
#' @param metric a name of FIGURES
#' @param view "rho05", "rho0" or "rhodiff"
surv_figure <- function(summary, metric, view) {
  spec <- FIGURES[[metric]]
  diff <- view == "rhodiff"

  # NA Pearson (a constant truth) has a NaN mean; leave it out rather than give
  # it a dodge slot
  keep <- filter(summary, is.finite(.data[[paste0("mean_", metric)]]))
  events_only <- paste0(view, "_", metric) %in% EVENTS_ONLY
  if (events_only) keep <- filter(keep, target != "Combined")

  fig <- point_range_plot(
    keep,
    metric,
    if (diff) paste0(spec$diff_lab, "\n", DIFF_LAB) else spec$y_lab,
    x = "censoring",
    colour = "model",
    shape = "estimand",
    shape_palette = scale_shape_manual(values = ESTIMAND_SHAPES),
    facet_rows = "target",
    facet_cols = "scenario",
    facet_scales = "free_y",
    hline = if (diff) 0 else spec$hline,
    dodge_width = 0.8,
    point_size = 1.5,
    ci_alpha = 0.7
  ) +
    labs(x = "Censoring", colour = "Model", shape = "Estimand") +
    theme(
      axis.text.x = element_text(angle = 0),
      legend.position = "bottom",
      legend.box = "vertical"
    ) +
    guides(colour = guide_legend(nrow = 2))

  save_fig(
    paste0("surv_", view, "_", metric, ".png"),
    fig_path,
    width = FIG_WIDTH,
    height = if (events_only) EVENTS_ONLY_HEIGHT else FIG_HEIGHT,
    plot = fig
  )
  fig
}

# kept in `figs` for viewing interactively, figs[[view]][[metric]]
figs <- list()
for (view in c(names(RHOS), "rhodiff")) {
  runs <- view_runs(view)
  if (nrow(runs) == 0) {
    message("no ", view, " metrics - figures skipped")
    next
  }
  summary <- summarise_view(runs)
  figs[[view]] <- lapply(setNames(names(FIGURES), names(FIGURES)),
                         surv_figure, summary = summary, view = view)
}
