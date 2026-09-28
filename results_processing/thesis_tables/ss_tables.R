##########
# title: LaTeX tables for the thesis chapter - sample size, both outcomes
##########
# One table for the reported scenarios 1-4 and one for the supplementary 5-10,
# per outcome, each metric as "mean (MCSE)". Labels and the summary come from
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
# - Rejection rates are 0/1 per run, so they get the binomial MCSE. Runs with an
#   NA p-value (constant tau, or a test not run - indep_po for the causal
#   forest and T-learners) are left out of the denominator, as in the figures.

library(here)
source(here("R", "figures.R"))
source(here("R", "tables.R"))

ALPHA <- 0.05

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

tests <- c("BLP_p", "BLP_p_os", "indep_cate", "indep_po")

# the table's columns, left to right. Drop a row here to drop a column.
table_cols <- function(o) {
  est <- "Estimation"
  # no "\\%" or other macros in a spanning header: add_header_above() strips
  # one level of backslash even with escape = FALSE. The caption gives alpha.
  rej <- "Rejection rate"
  # two-line headers (header2 "" for one line) keep the columns as narrow as
  # their cells, so the table fits a landscape A4 page
  tibble::tribble(
    ~stem,             ~header1,   ~header2,          ~group, ~digits,
    "ate_bias",        "Bias",     "",                est,    o$digits_bias,
    "rel_ate_bias",    "Rel. bias", "(\\%)",          est,    1,
    "mse",             "MSE",      "",                est,    o$digits_mse,
    "rmse",            "RMSE",     "",                est,    3,
    "corr",            "Pearson",  "",                est,    2,
    "spearman",        "Spearman", "",                est,    2,
    "sign_acc",        "Sign",     "accuracy",        est,    2,
    "rej_BLP_p",       "BLP",      "(2-sided)",       rej,    3,
    "rej_BLP_p_os",    "BLP",      "(1-sided, HC3)",  rej,    3,
    "rej_indep_cate",  "CATE",     "indep.",          rej,    3,
    "rej_indep_po",    "PO",       "indep.",          rej,    3
  ) %>%
    mutate(header = ifelse(
      header2 == "", header1,
      paste0("\\shortstack[r]{", header1, "\\\\", header2, "}")
    ))
}

# runs per cell, read off the data rather than the design, so failed runs show
runs_text <- function(metrics) {
  r <- range(count(metrics, scenario, n, model, name = "runs")$runs)
  if (r[1] == r[2]) r[1] else paste0(r[1], "--", r[2])
}

short_caption_text <- function(o, scenarios) {
  paste0("Sample-size study, ", o$dir, " outcome, scenarios ", scenarios)
}

caption_text <- function(o, scenarios, runs) {
  paste0(
    "Sample-size study, ", o$dir, " outcome, scenarios ", scenarios,
    ": mean (Monte Carlo SE) over ", runs, " runs per cell. ",
    "Bias and relative bias are for the ATE; rejection rates are at ",
    "$\\alpha = ", ALPHA, "$ (type I error in the null scenario, power ",
    "otherwise). --- : not defined (correlation under no heterogeneity; the PO ",
    "independence test for estimators that share another's pseudo-outcome)."
  )
}

for (o in outcomes) {
  metrics <- readRDS(file.path(dirname(path), "results", o$dir,
                               paste0(o$prefix, "_metrics.RDS")))

  metrics <- metrics %>%
    mutate(
      rel_ate_bias = 100 * rel_ate_bias,
      corr = if_else(scenario == 1, NA_real_, corr),
      spearman = if_else(scenario == 1, NA_real_, spearman),
      across(all_of(tests), ~ as.numeric(as.numeric(.x) < ALPHA),
             .names = "rej_{.col}")
    )

  cols <- table_cols(o)
  rej_stems <- grep("^rej_", cols$stem, value = TRUE)

  metrics_summary <- summarise_metrics(
    metrics,
    c("scenario", "n", "model"),
    cols = setNames(cols$stem, cols$stem),
    binomial = rej_stems,
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
          o, t$text, runs_text(filter(metrics, scenario %in% t$scenarios))
        ),
        label = paste0(o$prefix, "_ss_", nm)
      )
    out_file <- file.path(tab_path, paste0(o$prefix, "_ss_", nm, ".tex"))
    writeLines(tex, out_file)
    message("wrote ", out_file)
  }
}
