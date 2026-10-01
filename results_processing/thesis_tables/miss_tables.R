##########
# title: LaTeX tables for the thesis appendix - missing covariates, both outcomes
##########
# The full grid, one metric per table: rows are mechanism / scenario / handling
# method, the columns the five estimators, each cell "mean (MCSE)". That folds
# the summary's 560 rows (14 scenario x mechanism cells x 40 method x model
# cells) into 126 per table. The HTE tests get one table per test, with a True
# CATE column. Labels and the summary come from R/figures.R, the layout from
# R/tables.R (miss_latex_table()).
#
# Writes to ../results/thesis_tables/, for each outcome prefix (cts, bin):
#   <prefix>_miss_<metric>.tex      one per row of table_cols() below
#   <prefix>_miss_tests_<test>.tex  one per HTE test
# \input{} them into a document with
#   \usepackage{booktabs, longtable, pdflscape, array}
#
# Metric choices (missing/ADEMP.md, "Performance measures"):
# - Every estimation metric is on the complete units (_cu), the units every
#   handling method has. The all-unit versions compare complete case and IPW
#   (~350 units) with the rest (500), which they cannot.
# - Relative efficiency and the bias difference are against the complete-data
#   arm, where they are 1 and 0 by construction, so its rows are left out.
# - RMSE on the incomplete units is the error floor of scoring against tau at
#   covariates no method saw. Complete case and IPW have no incomplete units,
#   so their rows drop out.
# - RMSE against the truth given completeness is MNAR-tau's secondary truth,
#   defined on MNAR-tau rows only.
# - Correlation is undefined in scenario 1 (stored as 0, R/metrics.R), so it is
#   set to NA and prints as a dash.
# - A dash elsewhere is an estimator not fitted for that method: no DR-oracle
#   or DR-SuperLearner under multiple imputation, only the two grf estimators
#   under inbuilt missingness handling.

library(here)
source(here("R", "figures.R"))
source(here("R", "tables.R"))

ALPHA <- 0.05

# paths
path <- here()
res_root <- file.path(dirname(path), "results")
tab_path <- file.path(res_root, "thesis_tables")
dir.create(tab_path, showWarnings = FALSE, recursive = TRUE)

# the binary outcome's errors are on the risk-difference scale, an order of
# magnitude smaller, so bias gets an extra decimal
outcomes <- list(
  list(dir = "continuous", prefix = "cts", digits_bias = 3),
  list(dir = "binary", prefix = "bin", digits_bias = 4)
)

# one table per row. `what` opens the caption and is the list-of-tables entry,
# `note` follows the cell description. Drop a row here to drop a table.
table_cols <- function(o) {
  tibble::tribble(
    ~stem,                    ~digits,       ~drop_ref, ~what,                                                   ~note,
    "bias_cu",                o$digits_bias, FALSE,     "Bias, complete units",                                  "",
    "rel_ate_bias_cu",        1,             FALSE,     "Relative ATE bias (\\%), complete units",               "Bias of the complete units' ATE as a percentage of their true ATE.",
    "rmse_cu",                3,             FALSE,     "RMSE, complete units",                                  "",
    "corr_cu",                2,             FALSE,     "Pearson correlation with the true CATE, complete units", "Not defined in the null scenario.",
    "rel_efficiency_cu",      2,             TRUE,      "Relative efficiency, complete units",                   "MSE over the complete-data arm's MSE, same run and estimator: 1 is no loss.",
    "bias_diff_complete_cu",  o$digits_bias, TRUE,      "Bias difference from complete data, complete units",    "Bias minus the complete-data arm's bias, same run and estimator: the bias the missingness and its handling add.",
    "rmse_iu",                3,             FALSE,     "RMSE, incomplete units",                                "Complete case and IPW have no incomplete units.",
    "rmse_r_cu",              3,             FALSE,     "RMSE against the truth given completeness, complete units", "MNAR-$\\tau$ only: the true CATE plus the mean of $U$'s term among the complete units (the secondary truth)."
  )
}

# p-value column -> table title, in output order
tests <- c(
  BLP_p = "BLP (two-sided)",
  BLP_p_os = "BLP (one-sided, HC3)",
  indep_cate = "CATE independence",
  indep_po = "PO independence"
)

# p-values -> 0/1 rejection indicators, rej_<test>. A test absent from `df`
# (indep_po in the true-CATE file) is left absent.
add_rejections <- function(df) {
  df %>%
    mutate(across(any_of(names(tests)), ~ as.numeric(as.numeric(.x) < ALPHA),
                  .names = "rej_{.col}"))
}

# MECHANISM_LABELS spells MNAR-tau with a Unicode tau, which pdflatex cannot
# typeset
tex_mechanism <- function(df) {
  lv <- levels(df$mechanism)
  levels(df$mechanism)[lv == MECHANISM_LABELS[["MNAR-tau"]]] <- "MNAR-$\\tau$"
  df
}

group_cols <- c("scenario", "mechanism", "method")
run_cols <- c(group_cols, "model")

caption_text <- function(o, what, note, runs) {
  paste0(
    what, ", missing-covariate study, ", o$dir, " outcome: mean (Monte Carlo ",
    "SE) over ", runs, " runs per cell. ", note, if (note != "") " ",
    "Complete units: those with no amputed covariate. --- : estimator not ",
    "fitted for that handling method, or metric not defined."
  )
}

test_caption_text <- function(o, test, runs) {
  paste0(
    test, " test, missing-covariate study, ", o$dir, " outcome: rejection ",
    "rate at $\\alpha = ", ALPHA, "$ (Monte Carlo SE) over ", runs, " runs per ",
    "cell. Type I error in the null scenario, power otherwise. Runs with an ",
    "undefined p-value are left out, and handling methods with none at all. ",
    "True CATE: the test applied to the true CATE and nuisances, at the ",
    "handled covariates. --- : not run or not defined."
  )
}

for (o in outcomes) {
  res_path <- file.path(res_root, "missing", o$dir)
  metrics <- readRDS(file.path(res_path, paste0(o$prefix, "_miss_metrics.RDS")))
  runs <- runs_per_cell(metrics, run_cols)

  # --- estimation metrics -----------------------------------------------------
  cols <- table_cols(o)

  metrics_summary <- metrics %>%
    mutate(
      rel_ate_bias_cu = 100 * rel_ate_bias_cu,
      corr_cu = if_else(scenario == 1, NA_real_, corr_cu)
    ) %>%
    summarise_metrics(run_cols, cols = setNames(cols$stem, cols$stem),
                      count_na = character()) %>%
    apply_labels(MISS_SCENARIO_LABELS) %>%
    tex_mechanism()

  for (i in seq_len(nrow(cols))) {
    col <- cols[i, ]
    if (!paste0("mean_", col$stem) %in% names(metrics_summary)) {
      warning(o$prefix, "_miss_metrics.RDS has no ", col$stem, " - skipping")
      next
    }
    summ <- metrics_summary
    if (col$drop_ref) {
      summ <- filter(summ, method != METHOD_LABELS[["complete_data"]])
    }
    tex <- miss_latex_table(
      droplevels(summ), col$stem, col$digits,
      caption.short = paste0(col$what, ", missing covariates, ", o$dir, " outcome"),
      caption = caption_text(o, col$what, col$note, runs),
      label = paste0(o$prefix, "_miss_", col$stem)
    )
    out_file <- file.path(tab_path, paste0(o$prefix, "_miss_", col$stem, ".tex"))
    writeLines(tex, out_file)
    message("wrote ", out_file)
  }

  # --- HTE tests --------------------------------------------------------------
  true_cate_tests <- readRDS(
    file.path(res_path, paste0(o$prefix, "_miss_true_cate_tests.RDS"))
  ) %>%
    add_rejections()

  rej_cols <- paste0("rej_", names(tests))
  summarise_rejections <- function(df, group_cols) {
    summarise_metrics(df, group_cols,
                      cols = setNames(rej_cols, rej_cols)[rej_cols %in% names(df)],
                      binomial = rej_cols, count_na = character())
  }

  # apply_labels() puts the unlabelled "True CATE" after the estimators
  test_summary <- bind_rows(
    summarise_rejections(add_rejections(metrics), run_cols),
    summarise_rejections(true_cate_tests, group_cols) %>%
      mutate(model = "True CATE")
  ) %>%
    apply_labels(MISS_SCENARIO_LABELS) %>%
    tex_mechanism()

  for (test in names(tests)) {
    tex <- miss_latex_table(
      test_summary, paste0("rej_", test), 3,
      caption.short = paste0(tests[[test]], " test, missing covariates, ",
                             o$dir, " outcome"),
      caption = test_caption_text(o, tests[[test]], runs),
      label = paste0(o$prefix, "_miss_tests_", test)
    )
    out_file <- file.path(tab_path, paste0(o$prefix, "_miss_tests_", test, ".tex"))
    writeLines(tex, out_file)
    message("wrote ", out_file)
  }
}
