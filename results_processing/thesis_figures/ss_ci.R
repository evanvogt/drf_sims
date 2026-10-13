##########
# title: figures for the thesis chapter - sample size, CI studies
##########
# The correlated-covariate CI studies (sample_size/correlated/
# confidence_intervals/), per outcome, the CI_sf sweep: one figure per
# interval metric (marginal coverage, simultaneous coverage, mean interval
# length) at rho = 0, at rho = 0.5 and as the paired rho = 0.5 - rho = 0
# difference. The scenarios are across, the sample sizes down, the
# subsampling ratio on x and the models in colour. The tables are
# thesis_tables/ss_ci_tables.R. The data-driven CI_sf study (optimal_sf/)
# will be added here once its results are in.
#
# - Only the per-unit intervals, as the tables: the half-sample bootstrap band
#   of each model, and the causal forest's own variance-based interval
#   (causal_forest_inbuilt). The query-grid rows (*_grid) are left out.
# - The inbuilt interval does not depend on CI_sf (ss_ci_tables.R), so it is
#   drawn as a flat line across each panel, its 95% MCSE interval shaded,
#   from its rows at the first ratio.
# - Coverage shares one y scale across all panels; length is free per panel,
#   as the scenarios' CATE spread differs. facet_grid() can only free y per
#   row, so a length figure is one plot per scenario, put side by side.
# - Simultaneous coverage is 0/1 per run, so its levels get the binomial MCSE,
#   marginal coverage the general one. Differences are per run
#   (paired_rho_diff()), MCSE sd(diff) / sqrt(pairs).
# - An outcome whose metrics are missing, or do not yet cover every cell of
#   the design, is skipped with a message rather than drawn in part.
#
# Writes to ../results/thesis_figures/sample_size/:
#   {cts,bin}_corr_{rho0,rho05,rhodiff}_ci_{marg,simul,len}.png

library(here)
library(patchwork)
source(here("R", "figures.R"))

# paths
path <- here()
res_root <- file.path(dirname(path), "results", "correlated", "confidence_intervals")
fig_path <- file.path(dirname(path), "results", "thesis_figures", "sample_size")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

N_LEVELS <- c(500, 1000)
CI_SF <- round(seq(0.05, 0.5, 0.05), 2)
NOMINAL <- 0.95

# the study's design (*_corr_ci_config.R), to check the metrics are complete
DESIGN <- expand.grid(rho = c(0, 0.5), scenario = 1:4, n = N_LEVELS, CI_sf = CI_SF)

# second line of a difference figure's y labels
DIFF_LAB <- "ρ = 0.5 − ρ = 0"

# the one interval that does not depend on CI_sf
INBUILT <- "causal_forest_inbuilt"

# one figure each: `file` names the output, `binomial` marks a per-run 0/1
# indicator, `coverage` gets the nominal line and a fixed y scale
CI_METRICS <- tibble::tribble(
  ~stem,       ~col,                    ~file,   ~lab,                    ~coverage, ~binomial,
  "marg_cov",  "marginal_coverage",     "marg",  "Marginal coverage",     TRUE,      FALSE,
  "simul_cov", "simultaneous_coverage", "simul", "Simultaneous coverage", TRUE,      TRUE,
  "ci_len",    "mean_ci_length",        "len",   "Mean interval length",  FALSE,     FALSE
)

OUTCOMES <- list(
  list(dir = "continuous", prefix = "cts", scale_lab = ""),
  list(dir = "binary", prefix = "bin", scale_lab = "risk difference")
)

VIEWS <- list(
  rho0 = list(rho = 0, diff = FALSE),
  rho05 = list(rho = 0.5, diff = FALSE),
  rhodiff = list(diff = TRUE)
)

#' A metrics file, or NULL with a message if it is missing or incomplete
read_complete <- function(file, what) {
  if (!file.exists(file)) {
    message(what, ": ", file, " not found - skipped")
    return(NULL)
  }
  metrics <- readRDS(file)
  # CI_sf is stored as seq()'s doubles (0.15 is 0.15000000000000002), so it
  # is rounded here and in CI_SF, and the join is exact
  metrics <- mutate(metrics, CI_sf = round(CI_sf, 2))
  by <- c("rho", "scenario", "n", "CI_sf")
  gap <- anti_join(DESIGN, distinct(metrics, across(all_of(by))), by = by)
  if (nrow(gap) > 0) {
    message(what, ": ", nrow(gap), " of ", nrow(DESIGN), " design cells have ",
            "no metrics yet - skipped")
    return(NULL)
  }
  metrics
}

#' Mean and MCSE of each interval metric by scenario, n, CI_sf and model, at
#' one rho or as the paired difference
ci_summary <- function(metrics, view) {
  keys <- c("scenario", "n", "CI_sf", "model")
  metrics <- if (view$diff) {
    paired_rho_diff(metrics, CI_METRICS$col, c(keys, "run"))
  } else {
    filter(metrics, rho == view$rho)
  }
  summarise_metrics(
    metrics, keys,
    cols = setNames(CI_METRICS$col, CI_METRICS$stem),
    binomial = if (view$diff) character() else CI_METRICS$stem[CI_METRICS$binomial],
    count_na = character()
  ) %>%
    apply_labels(SS_SCENARIO_LABELS) %>%
    # n is a strip here, not the x axis, so it says what it is
    mutate(
      n = factor(n, levels = N_LEVELS, labels = paste("n =", N_LEVELS)),
      CI_sf = factor(CI_sf, levels = CI_SF)
    )
}

#' One metric by subsampling ratio, rows n, columns scenario
#'
#' @param summary output of ci_summary()
#' @param m one row of CI_METRICS
#' @param scale_lab the y label's last line, in brackets; "" for none
#' @param diff TRUE for a difference figure
ci_figure <- function(summary, m, scale_lab = "", diff = FALSE) {
  inbuilt_lab <- MODEL_LABELS[[INBUILT]]
  inbuilt <- filter(summary, model == inbuilt_lab)
  bands <- filter(summary, model != inbuilt_lab)

  # each model keeps the colour it has in the other sample-size figures (Safe,
  # by its place in MODEL_LABELS), and the inbuilt interval's shading takes its
  # line's colour. The limits fix the legend order: the inbuilt layers are
  # drawn first, so would otherwise put it first
  model_levels <- levels(summary$model)
  safe <- as.character(paletteer_d("rcartocolor::Safe", length(MODEL_LABELS)))
  pal_values <- setNames(safe[match(model_levels, MODEL_LABELS)], model_levels)

  lead <- if (diff) paste("Δ", tolower(m$lab)) else m$lab
  y_lab <- paste(
    c(lead, if (diff) DIFF_LAB, if (m$col == "mean_ci_length" && scale_lab != "")
      paste0("(", scale_lab, ")")),
    collapse = "\n"
  )
  hline <- if (diff) 0 else if (m$coverage) NOMINAL else NULL

  mean_col <- paste0("mean_", m$stem)
  mcse_col <- paste0("mcse_", m$stem)
  z <- qnorm(0.975)

  panel <- function(bands, inbuilt, facet_scales) {
    p <- point_range_plot(
      bands,
      m$stem,
      y_lab,
      x = "CI_sf",
      colour = "model",
      facet_rows = "n",
      facet_cols = "scenario",
      facet_scales = facet_scales,
      line = TRUE,
      ci_alpha = 0.7,
      hline = hline,
      palette = scale_colour_manual(values = pal_values, limits = model_levels)
    )
    # under the bands' points and lines
    p$layers <- c(
      list(
        geom_rect(
          data = inbuilt,
          aes(
            xmin = -Inf,
            xmax = Inf,
            ymin = .data[[mean_col]] - z * .data[[mcse_col]],
            ymax = .data[[mean_col]] + z * .data[[mcse_col]],
            fill = model
          ),
          alpha = 0.2,
          inherit.aes = FALSE
        ),
        geom_hline(
          data = inbuilt,
          aes(yintercept = .data[[mean_col]], colour = model),
          linewidth = 0.5
        )
      ),
      p$layers
    )
    p +
      scale_fill_manual(values = pal_values, limits = model_levels, guide = "none") +
      labs(x = "Subsampling ratio", colour = "Model")
  }

  fig <- if (m$coverage) {
    wrap_plots(panel(bands, inbuilt, "fixed"))
  } else {
    scenarios <- levels(droplevels(bands$scenario))
    plots <- lapply(seq_along(scenarios), function(i) {
      p <- panel(
        filter(bands, scenario == scenarios[i]),
        filter(inbuilt, scenario == scenarios[i]),
        "free_y"
      )
      # the n strips once, on the right
      if (i < length(scenarios)) {
        p <- p + theme(strip.text.y = element_blank(), strip.background.y = element_blank())
      }
      p
    })
    wrap_plots(plots, nrow = 1) +
      plot_layout(axis_titles = "collect")
  }

  fig +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom") &
    guides(colour = guide_legend(nrow = 2))
}

# ---- the CI_sf sweep ----------------------------------------------------------
# The figures are kept in `figs` for viewing interactively, named
# <prefix>_<rho0|rho05|rhodiff>_<marg|simul|len>.

figs <- list()
for (o in OUTCOMES) {
  metrics <- read_complete(
    file.path(res_root, o$dir, paste0(o$prefix, "_corr_ci_metrics.RDS")),
    o$prefix
  )
  if (is.null(metrics)) next

  metrics <- filter(metrics, !grepl("_grid$", model))

  # the inbuilt interval should be the same at every ratio; if a run's isn't,
  # the first ratio's still stands for it, but say so
  varies <- metrics %>%
    filter(model == INBUILT) %>%
    group_by(rho, scenario, n, run) %>%
    summarise(
      across(all_of(CI_METRICS$col), ~ n_distinct(round(.x, 10)) > 1),
      .groups = "drop"
    ) %>%
    filter(if_any(all_of(CI_METRICS$col)))
  if (nrow(varies) > 0) {
    warning(o$prefix, ": the inbuilt interval varies with CI_sf in ",
            nrow(varies), " runs - drawn from CI_sf = ", CI_SF[1])
  }
  metrics <- filter(metrics, model != INBUILT | CI_sf == CI_SF[1])

  for (v in names(VIEWS)) {
    summary <- ci_summary(metrics, VIEWS[[v]])
    for (i in seq_len(nrow(CI_METRICS))) {
      m <- CI_METRICS[i, ]
      name <- paste(o$prefix, v, m$file, sep = "_")
      figs[[name]] <- ci_figure(summary, m, o$scale_lab, VIEWS[[v]]$diff)
      save_fig(
        paste0(o$prefix, "_corr_", v, "_ci_", m$file, ".png"),
        fig_path,
        width = 21,
        height = 13,
        plot = figs[[name]]
      )
    }
  }
}
