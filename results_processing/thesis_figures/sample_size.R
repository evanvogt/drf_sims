##########
# title: figures for the thesis chapter - sample size, correlated covariates
##########
# The sample-size chapter's figures, from the correlated-covariate studies
# (sample_size/correlated/). Per outcome, for the main text's scenarios 1-4
# and the appendix's 5-10 (*_supp_*), each at rho = 0.5 and as the paired
# rho = 0.5 - rho = 0 difference (*_rhodiff_*):
#   - CATE bias on the top row, RMSE on the bottom, the scenarios across;
#   - each HTE test's rejection rate at HTE_ALPHA, one row per test, with the
#     same tests run on the true values as a reference series.
# The tables are thesis_tables/ss_tables.R and ss_test_tables.R.
#
# The differences: run r at rho = 0 and at rho = 0.5 shares its random draws
# (sample_size/correlated/README.md), so each metric is differenced per run
# (paired_rho_diff()) and its MCSE is sd(diff) / sqrt(pairs). The truth moves
# with rho too - bW, the ATE, SD(tau) and cor(m0, tau)
# (sample_size/correlated/corr_truth_summary.R) - so a difference is the
# models' response to the whole change in DGM, not to the correlation alone.
#
# Labels, palette, summaries, HTE_ALPHA and figure sizing come from
# R/figures.R. This script carries only the paths and this study's filters.

library(here)
library(patchwork)
source(here("R", "figures.R"))

# paths
path <- here()
res_path <- file.path(dirname(path), "results", "correlated")
fig_path <- file.path(dirname(path), "results", "thesis_figures", "sample_size")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

RHO <- 0.5
N_LEVELS <- c(100, 250, 500, 1000)

# second line of a difference figure's y labels
DIFF_LAB <- "ρ = 0.5 − ρ = 0"

# the appendix's scenarios: six panels across rather than four, so wider
SUPP_SCENARIOS <- 5:10
SUPP_WIDTH <- 29

#' One of a study's per-run files, both rhos, the given scenarios
read_corr <- function(outcome, prefix, stem, scenarios) {
  readRDS(file.path(
    res_path,
    outcome,
    paste0(prefix, "_corr_", stem, ".RDS")
  )) %>%
    filter(scenario %in% scenarios)
}

#' Bias (top) and RMSE (bottom) by sample size, one panel per scenario
#'
#' @param outcome results folder, "continuous" or "binary"
#' @param prefix file prefix, "cts" or "bin"
#' @param scale_lab the y labels' last line, in brackets, e.g. "risk
#'   difference"; "" for none
#' @param scenarios,labels the scenarios to plot, and their strip labels
#' @param suffix inserted into the file name, e.g. "_supp"
#' @param width figure width in cm, as save_fig()
#' @param diff FALSE for the rho = RHO levels, TRUE for the paired
#'   rho = 0.5 - rho = 0 differences
bias_rmse_figure <- function(
  outcome,
  prefix,
  scale_lab = "",
  scenarios = 1:4,
  labels = SS_SCENARIO_LABELS,
  suffix = "",
  width = 21,
  diff = FALSE
) {
  metrics <- read_corr(outcome, prefix, "metrics", scenarios)
  metrics <- if (diff) {
    paired_rho_diff(metrics, c("bias", "rmse"))
  } else {
    filter(metrics, rho == RHO)
  }
  metrics <- metrics %>%
    apply_labels(labels) %>%
    mutate(n = factor(n, levels = N_LEVELS))

  # rmse isn't in summarise_metrics()'s default cols
  metrics_summary <- summarise_metrics(
    metrics,
    c("scenario", "n", "model"),
    cols = c(bias = "bias", rmse = "rmse")
  )

  # one short line each, or the two rows' labels run into each other
  y_lab <- function(what) {
    # "bias" or "RMSE", capitalised to start a levels label
    lead <- if (diff) paste("Δ", what) else sub("^b", "B", what)
    first <- paste(lead, "of the CATE")
    paste(
      c(
        first,
        if (diff) DIFF_LAB,
        if (scale_lab != "") paste0("(", scale_lab, ")")
      ),
      collapse = "\n"
    )
  }

  bias_row <- point_range_plot(
    metrics_summary,
    "bias",
    y_lab("bias"),
    x = "n",
    colour = "model",
    facet_rows = NULL,
    facet_cols = "scenario",
    line = TRUE,
    ci_alpha = 0.7,
    hline = 0,
    blank_x = TRUE
  ) +
    labs(x = NULL)

  # scenario names are on the top row's strips already. No change is the
  # RMSE difference's reference; the levels have none
  rmse_row <- point_range_plot(
    metrics_summary,
    "rmse",
    y_lab("RMSE"),
    x = "n",
    colour = "model",
    facet_rows = NULL,
    facet_cols = "scenario",
    line = TRUE,
    ci_alpha = 0.7,
    hline = if (diff) 0 else NULL
  ) +
    labs(x = "Sample size") +
    theme(strip.text.x = element_blank(), strip.background.x = element_blank())

  fig <- (bias_row / rmse_row) +
    plot_layout(guides = "collect") &
    labs(colour = "Model") &
    theme(legend.position = "bottom") &
    guides(colour = guide_legend(nrow = 2))
  save_fig(
    paste0(prefix, "_corr", suffix, if (diff) "_rhodiff", "_bias_rmse.png"),
    fig_path,
    width = width,
    plot = fig
  )
  fig
}

# the power each HTE test should reach where there is HTE (a conventional
# target - nothing in the DGM aims for it)
POWER_TARGET <- 0.8

# the HTE tests, in facet-row order
HTE_TESTS <- c(
  BLP_p_os = "BLP (one-sided, HC3)",
  BLP_p = "BLP (two-sided)",
  indep_po = "PO independence",
  indep_cate = "CATE independence"
)

TRUE_LAB <- "True values"

#' Rejection rate at HTE_ALPHA by sample size, one row per test, one column
#' per scenario
#'
#' Levels: dashed reference line at HTE_ALPHA in the null scenario (nominal
#' size) and at POWER_TARGET in the others. Differences: the per-run
#' difference in the 0/1 rejection, so each point is the change in rejection
#' rate, with a dashed line at 0 throughout.
#'
#' The "True values" series is the BLP and CATE independence tests run on the
#' true CATE and true nuisances (<prefix>_corr_true_cate_tests.RDS, see
#' true_cate_test_row() in R/cate_models.R). The true CATE is constant in the
#' null scenario, so those tests have no true series there. There is no
#' true-value PO test: DR-oracle's own indep_po already tests the true
#' pseudo-outcome, so it stands in as the PO row's true series rather than as a
#' model of its own. Causal forest and the T-learners have no indep_po (see
#' hte_test_metrics() in R/metrics.R), so they draw nothing in that row. In a
#' difference figure the true series is the reference: the DGM's own power
#' moves with rho.
#'
#' @param outcome results folder, "continuous" or "binary"
#' @param prefix file prefix, "cts" or "bin"
#' @param scenarios,labels,suffix,width,diff as bias_rmse_figure()
rejection_figure <- function(
  outcome,
  prefix,
  scenarios = 1:4,
  labels = SS_SCENARIO_LABELS,
  suffix = "",
  width = 21,
  diff = FALSE
) {
  metrics <- read_corr(outcome, prefix, "metrics", scenarios)
  true_tests <- read_corr(outcome, prefix, "true_cate_tests", scenarios)

  est <- metrics %>%
    transmute(
      rho,
      run,
      scenario,
      n,
      model = as.character(model),
      across(all_of(names(HTE_TESTS)))
    ) %>%
    pivot_longer(
      all_of(names(HTE_TESTS)),
      names_to = "test",
      values_to = "p"
    ) %>%
    filter(!(model == "dr_oracle" & test == "indep_po"))

  true_cols <- c("BLP_p", "BLP_p_os", "indep_cate")
  truth <- bind_rows(
    true_tests %>%
      select(rho, run, scenario, n, all_of(true_cols)) %>%
      pivot_longer(all_of(true_cols), names_to = "test", values_to = "p"),
    metrics %>%
      filter(model == "dr_oracle") %>%
      transmute(rho, run, scenario, n, test = "indep_po", p = indep_po)
  ) %>%
    mutate(model = TRUE_LAB)

  rejections <- bind_rows(est, truth) %>%
    mutate(rej = as.numeric(p < HTE_ALPHA))
  rejections <- if (diff) {
    paired_rho_diff(
      rejections,
      "rej",
      keys = c("scenario", "n", "model", "test", "run")
    )
  } else {
    filter(rejections, rho == RHO)
  }

  # NA p-values (a constant CATE, or a test with no value for that model) are
  # left out, as summarise_metrics()'s na.rm would anyway; in a difference, NA
  # at either rho drops the pair
  rejections <- rejections %>%
    filter(!is.na(rej)) %>%
    mutate(test = factor(test, levels = names(HTE_TESTS), labels = HTE_TESTS)) %>%
    apply_labels(labels) %>%
    mutate(n = factor(n, levels = N_LEVELS))

  # a rate is a per-run 0/1 indicator, so the binomial MCSE; a difference of
  # two is not, and the general MCSE on it is the paired SE
  rej_summary <- summarise_metrics(
    rejections,
    c("test", "scenario", "n", "model"),
    cols = c(rej = "rej"),
    binomial = if (diff) character() else "rej",
    count_na = character()
  )

  # the models keep their colours from the bias/RMSE figure (Safe, in level
  # order); the true series, label_factor()'s last level, is black
  model_levels <- levels(rej_summary$model)
  pal <- scale_colour_manual(
    values = setNames(
      c(
        as.character(paletteer_d(
          "rcartocolor::Safe",
          length(model_levels) - 1
        )),
        "black"
      ),
      model_levels
    )
  )

  # levels: nominal size where there is no HTE, the power target everywhere
  # else (scenarios 5-10 all have HTE). Differences: no change
  scenario_levels <- levels(rej_summary$scenario)
  ref_lines <- tibble(
    scenario = factor(scenario_levels, levels = scenario_levels),
    yintercept = if (diff) {
      0
    } else {
      if_else(
        scenario_levels == SS_SCENARIO_LABELS[["1"]],
        HTE_ALPHA,
        POWER_TARGET
      )
    }
  )

  y_lab <- paste("Rejection rate at", HTE_ALPHA)
  if (diff) y_lab <- paste0("Δ ", tolower(y_lab), "\n", DIFF_LAB)

  fig <- point_range_plot(
    rej_summary,
    "rej",
    y_lab,
    x = "n",
    colour = "model",
    facet_rows = "test",
    facet_cols = "scenario",
    facet_scales = "fixed",
    line = TRUE,
    ci_alpha = 0.7,
    hline = ref_lines,
    palette = pal
  ) +
    labs(x = "Sample size", colour = "Model") +
    theme(legend.position = "bottom") +
    guides(colour = guide_legend(nrow = 2))
  save_fig(
    paste0(prefix, "_corr", suffix, if (diff) "_rhodiff", "_rejection.png"),
    fig_path,
    width = width,
    height = 18,
    plot = fig
  )
  fig
}

# ---- all figures -------------------------------------------------------------
# Per outcome, scenario set and view. The appendix's set is drawn only once an
# outcome's metrics have scenarios 5-10 (they were added to the study after 1-4
# had run). The figures are kept in `figs` for viewing interactively, named
# <prefix>_<set>_<rho05|rhodiff>.

OUTCOMES <- list(
  list(outcome = "continuous", prefix = "cts", scale_lab = ""),
  list(outcome = "binary", prefix = "bin", scale_lab = "risk difference")
)

FIG_SETS <- list(
  main = list(
    scenarios = 1:4,
    labels = SS_SCENARIO_LABELS,
    suffix = "",
    width = 21
  ),
  supp = list(
    scenarios = SUPP_SCENARIOS,
    labels = SS_SUPP_SCENARIO_PLOT_LABELS,
    suffix = "_supp",
    width = SUPP_WIDTH
  )
)

figs <- list()
for (o in OUTCOMES) {
  for (set in names(FIG_SETS)) {
    s <- FIG_SETS[[set]]
    has_set <- nrow(read_corr(o$outcome, o$prefix, "metrics", s$scenarios)) > 0
    if (!has_set) {
      message(o$outcome, ": no ", set, " scenario metrics yet - figures skipped")
      next
    }
    for (is_diff in c(FALSE, TRUE)) {
      figs[[paste(o$prefix, set, if (is_diff) "rhodiff" else "rho05", sep = "_")]] <-
        list(
          bias_rmse = bias_rmse_figure(
            o$outcome,
            o$prefix,
            scale_lab = o$scale_lab,
            scenarios = s$scenarios,
            labels = s$labels,
            suffix = s$suffix,
            width = s$width,
            diff = is_diff
          ),
          rejection = rejection_figure(
            o$outcome,
            o$prefix,
            scenarios = s$scenarios,
            labels = s$labels,
            suffix = s$suffix,
            width = s$width,
            diff = is_diff
          )
        )
    }
  }
}
