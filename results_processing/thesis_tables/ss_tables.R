##########
# title: LaTeX tables for the thesis chapter - sample size, both outcomes
##########
# One table for the reported scenarios 1-4 and one for the supplementary 5-10,
# per outcome, each estimation metric as "mean (MCSE)". The HTE tests have
# their own tables, ss_test_tables.R. Labels and the summary come from
# R/figures.R, the table layout from R/tables.R.
#
# Writes to ../results/thesis_tables/:
#   cts_ss_main.tex, cts_ss_supp.tex, bin_ss_main.tex, bin_ss_supp.tex
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
outcomes <- list(
  list(dir = "continuous", prefix = "cts", digits_bias = 3, digits_mse = 3),
  list(dir = "binary", prefix = "bin", digits_bias = 4, digits_mse = 4)
)

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
  paste0("CATE estimation, ", o$dir, " outcome, scenarios ", scenarios)
}

caption_text <- function(o, scenarios, runs) {
  paste0(
    "CATE estimation, sample-size study, ", o$dir, " outcome, scenarios ",
    scenarios, ": mean (Monte Carlo SE) over ", runs, " runs per cell. ",
    "Bias and relative bias are for the ATE. --- : not defined (correlation ",
    "under no heterogeneity)."
  )
}

for (o in outcomes) {
  metrics <- readRDS(file.path(dirname(path), "results", o$dir,
                               paste0(o$prefix, "_metrics.RDS")))

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

  for (nm in names(tables)) {
    t <- tables[[nm]]
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
        label = paste0(o$prefix, "_ss_", nm),
        # portrait fits the 7 columns within 16cm (A4, 2.5cm margins) with
        # 2pt column padding; the 4-decimal binary table needs it. Set
        # landscape = TRUE instead if the thesis margins are wider
        landscape = FALSE,
        tabcolsep = "2pt"
      )
    out_file <- file.path(tab_path, paste0(o$prefix, "_ss_", nm, ".tex"))
    writeLines(tex, out_file)
    message("wrote ", out_file)
  }
}
