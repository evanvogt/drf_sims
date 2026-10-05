##########
# title: LaTeX tables for the thesis chapter - sample size, both outcomes
##########
# The correlated-covariate studies (sample_size/correlated/), per outcome: the
# rho = 0.5 levels, and the paired rho = 0.5 - rho = 0 differences in place of
# the rho = 0 levels. Each as one table for the reported scenarios 1-4 and one
# for the supplementary 5-10, each estimation metric as "mean (MCSE)". The
# independent-covariate studies (sample_size/continuous/, binary/) are no
# longer tabulated. The HTE tests have their own tables, ss_test_tables.R.
# Labels and the summary come from R/figures.R, the table layout from
# R/tables.R.
#
# Writes to ../results/thesis_tables/:
#   {cts,bin}_corr_rho05_ss_main.tex, {cts,bin}_corr_rho05_ss_supp.tex,
#   {cts,bin}_corr_rhodiff_ss_main.tex, {cts,bin}_corr_rhodiff_ss_supp.tex
# \input{} them into a document with
#   \usepackage{booktabs, longtable, pdflscape, array}
#
# The differences: run r at rho = 0 and at rho = 0.5 shares its random draws
# (sample_size/correlated/README.md, "Seeding and pairing"), so each metric is
# differenced per run (paired_rho_diff() in R/figures.R) and its MCSE is
# sd(diff) / sqrt(pairs), as in correlated/*/*_corr_results.qmd's
# paired_diff(). Runs missing at either rho drop out of the pairs. Pearson
# and sign accuracy get a third decimal there: their differences are
# hundredths.
#
# Metric choices:
# - Bias is the CATE bias, `bias` (mean of est - true over units), as in the
#   figures. `ate_bias` is the same number unless a run has NA estimates, but
#   it is not an ATE estimate: that would need AIPW with the CATE models'
#   propensity scores, which not every model estimates.
# - No relative bias. The true ATE is recalibrated per n for 80% power and
#   shrinks with n, so the bias rows are not on a common scale across n.
# - MSE, not RMSE, and Pearson, not Spearman: the missing-data chapter's
#   tables report Pearson too (miss_tables.R), and in scenarios 5 and 7 the
#   true CATE is mostly ties, which muddies Spearman's ranks.
# - Pearson is undefined in scenario 1 (stored as 0, R/metrics.R), so it is
#   set to NA and prints as a dash.

library(here)
source(here("R", "figures.R"))
source(here("R", "tables.R"))

# paths
path <- here()
tab_path <- file.path(dirname(path), "results", "thesis_tables")
dir.create(tab_path, showWarnings = FALSE, recursive = TRUE)

# the binary outcome's errors are on the risk-difference scale, an order of
# magnitude smaller, so bias and MSE get an extra decimal. `supp_note` is
# appended to the 5-10 table's caption: binary scenario 9 runs on a smaller
# scale than binary/'s (RD_SCALE_CORR in R/dgm_scenarios.R)
outcomes <- list(
  list(dir = "continuous", prefix = "cts", label = "continuous outcome",
       digits_bias = 3, digits_mse = 3, supp_note = ""),
  list(dir = "binary", prefix = "bin", label = "binary outcome",
       digits_bias = 4, digits_mse = 4,
       supp_note = paste0(
         " In scenario 9 the heterogeneity is scaled by 0.139 rather than the",
         " independent-covariate study's 0.164, so that every treated risk",
         " stays within $[0.01, 0.99]$ at $\\rho = 0.5$."
       ))
)

# one entry per outcome and table set, the rho = 0.5 levels and the paired
# differences (CORR_RHOS in R/dgm_scenarios.R): `out` names the output files
# and LaTeX labels, `label` is the caption's outcome text
outcomes <- unlist(lapply(outcomes, function(o) list(
  modifyList(o, list(
    out = paste0(o$prefix, "_corr_rho05"),
    label = paste0(o$label, ", correlated covariates ($\\rho = 0.5$)"),
    diff = FALSE
  )),
  modifyList(o, list(
    out = paste0(o$prefix, "_corr_rhodiff"),
    label = paste0(o$label, ", correlated covariates, $\\rho = 0.5$ minus ",
                   "$\\rho = 0$"),
    diff = TRUE
  ))
)), recursive = FALSE)

# the table's columns, left to right. Drop a row here to drop a column.
table_cols <- function(o) {
  digits_prop <- if (o$diff) 3 else 2
  # two-line headers (header2 "" for one line) keep the columns as narrow as
  # their cells
  tibble::tribble(
    ~stem,             ~header1,    ~header2,    ~digits,
    "bias",            "Bias",      "",          o$digits_bias,
    "mse",             "MSE",       "",          o$digits_mse,
    "corr",            "Pearson",   "",          digits_prop,
    "sign_acc",        "Sign",      "accuracy",  digits_prop
  ) %>%
    mutate(header = ifelse(
      header2 == "", header1,
      paste0("\\shortstack[r]{", header1, "\\\\", header2, "}")
    ))
}

short_caption_text <- function(o, scenarios) {
  paste0("CATE estimation, ", o$label, ", scenarios ", scenarios)
}

caption_text <- function(o, scenarios, runs, note = "") {
  cells <- if (o$diff) {
    paste0(
      "per-run difference, mean (Monte Carlo SE) over ", runs, " run pairs per cell. A run shares its random draws ",
      "at the two values of $\\rho$, so the SE is that of the paired difference. "
    )
  } else {
    paste0("mean (Monte Carlo SE) over ", runs, " runs per cell. ")
  }
  paste0(
    "CATE estimation, sample-size study, ", o$label, ", scenarios ",
    scenarios, ": ", cells,
    "Bias is the CATE bias, the mean of the estimated minus true CATE over ",
    "units. --- : not defined (correlation ",
    "under no heterogeneity).", note
  )
}

# paired_rho_diff() is R/figures.R's, shared with ss_test_tables.R and the
# figures
for (o in outcomes) {
  metrics <- readRDS(file.path(dirname(path), "results", "correlated", o$dir,
                               paste0(o$prefix, "_corr_metrics.RDS")))

  metrics <- metrics %>%
    mutate(corr = if_else(scenario == 1, NA_real_, corr))

  cols <- table_cols(o)

  metrics <- if (o$diff) {
    paired_rho_diff(metrics, cols$stem)
  } else {
    filter(metrics, rho == 0.5)
  }

  metrics_summary <- summarise_metrics(
    metrics,
    c("scenario", "n", "model"),
    cols = setNames(cols$stem, cols$stem),
    count_na = character()
  )

  tables <- list(
    main = list(scenarios = 1:4, labels = SS_SCENARIO_LABELS,
                text = "1--4", note = ""),
    supp = list(scenarios = 5:10, labels = SS_SUPP_SCENARIO_LABELS,
                text = "5--10", note = o$supp_note)
  )

  for (nm in names(tables)) {
    t <- tables[[nm]]
    # 5-10 were added to the studies after 1-4 had run
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
          o, t$text, runs_per_cell(filter(metrics, scenario %in% t$scenarios)),
          t$note
        ),
        label = paste0(o$out, "_ss_", nm),
        # the 4 columns fit portrait within the thesis's 17cm (A4, 2cm side
        # margins) at 10pt with LaTeX's own column padding, the 4-decimal
        # binary differences included
        landscape = FALSE,
        font_size = 10,
        tabcolsep = NULL
      )
    out_file <- file.path(tab_path, paste0(o$out, "_ss_", nm, ".tex"))
    writeLines(tex, out_file)
    message("wrote ", out_file)
  }
}
