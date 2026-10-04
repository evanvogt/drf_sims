##########
# title: LaTeX tables for the thesis chapter - sample size, both outcomes
##########
# One table for the reported scenarios 1-4 and one for the supplementary 5-10,
# per outcome, each estimation metric as "mean (MCSE)". The correlated-covariate
# studies (sample_size/correlated/) get the supplementary table only, once per
# rho: their scenarios 1-4 are the main text's figures
# (thesis_figures/sample_size.R). The HTE tests have their own tables,
# ss_test_tables.R. Labels and the summary come from R/figures.R, the table
# layout from R/tables.R.
#
# Writes to ../results/thesis_tables/:
#   cts_ss_main.tex, cts_ss_supp.tex, bin_ss_main.tex, bin_ss_supp.tex,
#   {cts,bin}_corr_rho0_ss_supp.tex, {cts,bin}_corr_rho05_ss_supp.tex
# \input{} them into a document with
#   \usepackage{booktabs, longtable, pdflscape, array}
#
# Metric choices:
# - Bias is `ate_bias`. `bias` (mean of est - true over units) is the same
#   number unless a run has NA estimates, so it is not repeated.
# - Relative bias is the ATE's, in %. The true ATE is recalibrated per n for
#   80% power and shrinks with n, so this is what makes the n rows comparable.
#   The per-unit `rel_bias_cate` is left out: the true CATE crosses or
#   approaches zero in several scenarios, and those units dominate its mean.
# - Correlations are undefined in scenario 1 (stored as 0, R/metrics.R), so
#   they are set to NA and print as a dash.

library(here)
source(here("R", "figures.R"))
source(here("R", "tables.R"))

# paths
path <- here()
tab_path <- file.path(dirname(path), "results", "thesis_tables")
dir.create(tab_path, showWarnings = FALSE, recursive = TRUE)

# the binary outcome's errors are on the risk-difference scale, an order of
# magnitude smaller, so bias and MSE get an extra decimal
#
# Per entry: `dir` and `prefix` locate the metrics; `out` names the output
# files and LaTeX labels; `label` is the caption's outcome text; `rho` (the
# correlated studies only) is the rho kept; `tables` is which of main (1-4)
# and supp (5-10) to write; `note` is appended to the caption.
outcomes <- list(
  list(dir = "continuous", prefix = "cts", out = "cts",
       label = "continuous outcome", tables = c("main", "supp"), note = "",
       digits_bias = 3, digits_mse = 3),
  list(dir = "binary", prefix = "bin", out = "bin",
       label = "binary outcome", tables = c("main", "supp"), note = "",
       digits_bias = 4, digits_mse = 4)
)

# the correlated studies, at each rho (CORR_RHOS in R/dgm_scenarios.R).
# Binary scenario 9 runs on a smaller scale than binary/'s (RD_SCALE_CORR)
BIN_CORR_NOTE <- paste0(
  " In scenario 9 the heterogeneity is scaled by 0.139 rather than the",
  " independent-covariate study's 0.164, so that every treated risk stays",
  " within $[0.01, 0.99]$ at $\\rho = 0.5$."
)
corr_outcomes <- lapply(outcomes, function(o) {
  lapply(c(0, 0.5), function(r) modifyList(o, list(
    dir = file.path("correlated", o$dir),
    prefix = paste0(o$prefix, "_corr"),
    out = paste0(o$prefix, "_corr_rho", sub(".", "", r, fixed = TRUE)),
    label = paste0(o$label, ", correlated covariates ($\\rho = ", r, "$)"),
    rho = r,
    tables = "supp",
    note = if (o$dir == "binary") BIN_CORR_NOTE else ""
  )))
})
outcomes <- c(outcomes, unlist(corr_outcomes, recursive = FALSE))

# the table's columns, left to right. Drop a row here to drop a column.
table_cols <- function(o) {
  # two-line headers (header2 "" for one line) keep the columns as narrow as
  # their cells
  tibble::tribble(
    ~stem,             ~header1,    ~header2,    ~digits,
    "ate_bias",        "Bias",      "",          o$digits_bias,
    "rel_ate_bias",    "Rel. bias", "(\\%)",     1,
    "mse",             "MSE",       "",          o$digits_mse,
    "rmse",            "RMSE",      "",          3,
    "corr",            "Pearson",   "",          2,
    "spearman",        "Spearman",  "",          2,
    "sign_acc",        "Sign",      "accuracy",  2
  ) %>%
    mutate(header = ifelse(
      header2 == "", header1,
      paste0("\\shortstack[r]{", header1, "\\\\", header2, "}")
    ))
}

short_caption_text <- function(o, scenarios) {
  paste0("CATE estimation, ", o$label, ", scenarios ", scenarios)
}

caption_text <- function(o, scenarios, runs) {
  paste0(
    "CATE estimation, sample-size study, ", o$label, ", scenarios ",
    scenarios, ": mean (Monte Carlo SE) over ", runs, " runs per cell. ",
    "Bias and relative bias are for the ATE. --- : not defined (correlation ",
    "under no heterogeneity).", o$note
  )
}

for (o in outcomes) {
  metrics <- readRDS(file.path(dirname(path), "results", o$dir,
                               paste0(o$prefix, "_metrics.RDS")))
  if (!is.null(o$rho)) metrics <- filter(metrics, rho == o$rho)

  metrics <- metrics %>%
    mutate(
      rel_ate_bias = 100 * rel_ate_bias,
      corr = if_else(scenario == 1, NA_real_, corr),
      spearman = if_else(scenario == 1, NA_real_, spearman)
    )

  cols <- table_cols(o)

  metrics_summary <- summarise_metrics(
    metrics,
    c("scenario", "n", "model"),
    cols = setNames(cols$stem, cols$stem),
    count_na = character()
  )

  tables <- list(
    main = list(scenarios = 1:4, labels = SS_SCENARIO_LABELS,
                text = "1--4"),
    supp = list(scenarios = 5:10, labels = SS_SUPP_SCENARIO_LABELS,
                text = "5--10")
  )

  for (nm in o$tables) {
    t <- tables[[nm]]
    # the correlated studies' 5-10 were added after their 1-4 had run
    if (!any(metrics$scenario %in% t$scenarios)) {
      message(o$out, ": no scenario ", t$text, " metrics yet - ", nm, " table skipped")
      next
    }
    tex <- metrics_summary %>%
      filter(scenario %in% t$scenarios) %>%
      apply_labels(t$labels) %>%
      mutate(n = factor(n, levels = c(100, 250, 500, 1000))) %>%
      ss_latex_table(
        cols,
        caption.short = short_caption_text(o, t$text),
        caption = caption_text(
          o, t$text, runs_per_cell(filter(metrics, scenario %in% t$scenarios))
        ),
        label = paste0(o$out, "_ss_", nm),
        # portrait fits the 7 columns within 16cm (A4, 2.5cm margins) with
        # 2pt column padding; the 4-decimal binary table needs it. Set
        # landscape = TRUE instead if the thesis margins are wider
        landscape = FALSE,
        tabcolsep = "2pt"
      )
    out_file <- file.path(tab_path, paste0(o$out, "_ss_", nm, ".tex"))
    writeLines(tex, out_file)
    message("wrote ", out_file)
  }
}
