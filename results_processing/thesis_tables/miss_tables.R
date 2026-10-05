##########
# title: LaTeX tables for the thesis appendix - missing covariates, both outcomes
##########
# The main scenarios (1-4), one estimation table per estimator: rows are
# mechanism / scenario / handling method, the columns relative efficiency and
# CATE bias difference, each on all units, the complete units and the
# incomplete units, each cell "mean (MCSE)". The HTE tests get one table per
# test, the estimators across the columns plus a True CATE column. Labels and
# the summary come from R/figures.R, the layout from R/tables.R
# (miss_metrics_table(), miss_latex_table()).
#
# Writes to ../results/thesis_tables/, for each outcome prefix (cts, bin):
#   <prefix>_miss_<model>.tex       one per estimator, e.g. cts_miss_dr_random_forest.tex
#   <prefix>_miss_tests_<test>.tex  one per HTE test
# \input{} them into a document with
#   \usepackage{booktabs, longtable, pdflscape, array}
#
# Metric choices (missing/ADEMP.md, "Performance measures"):
# - Relative efficiency (MSE over the complete-data arm's MSE) and CATE bias
#   difference (CATE bias minus the complete-data arm's), same run and
#   estimator, as in the chapter's figures (thesis_figures/missing.R). The
#   complete-data arm is 1 and 0 by construction, so its rows are left out.
# - Each on three sets of units (cate_metrics_split() in R/metrics.R): all
#   units is what an analyst gets, the complete / incomplete split says where
#   the loss is. The complete-data arm is split by the amputation its run would
#   have had, so each set is compared with the same units in the reference.
#   Complete case and IPW drop the incomplete units, so they have only the
#   complete-unit columns: their all-unit versions would compare ~350 units
#   with 500, and are NA in the metrics file.
# - Bias is the CATE bias, mean of est - true over units, not an ATE bias:
#   that would need AIPW with the CATE models' propensity scores, which not
#   every model estimates.
# - No relative bias: the bias difference is the comparison instead.
#   `rel_bias_cate_cu` averages (est - true) / true per unit and blows up
#   wherever the true CATE is near 0, and `rel_ate_bias_cu` is relative to the
#   ATE.
# - No Pearson correlation. Split by units it would mislead: the complete
#   units' narrower range of true CATE lowers it under MAR even in the
#   complete-data arm.
# - The truth is the primary one, tau(X), under MNAR-tau too.
# - A handling method an estimator is not fitted for has no row in that
#   estimator's table: no DR-oracle or DR-SuperLearner under multiple
#   imputation, only the two grf estimators under inbuilt missingness handling.

library(here)
source(here("R", "figures.R"))
source(here("R", "tables.R"))

ALPHA <- 0.05

# the main scenarios: the main study's 1-4 (TE_MISS in R/dgm_scenarios.R), so
# they take the sample-size chapter's labels, as the chapter's figures do
SCENARIOS <- 1:4
REF_METHOD <- "complete_data"
# what a run shares with its complete-data reference (ref_by in
# missing/*/*_miss_metrics.R)
REF_BY <- c("scenario", "n", "type", "prop", "mechanism", "run", "model")

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

# the estimation tables' columns, left to right, under a spanning header per
# metric. Drop a row here to drop a column.
table_cols <- function(o) {
  tibble::tribble(
    ~stem,                   ~group,                 ~header,      ~digits,
    "rel_efficiency",        "Relative efficiency",  "All",        2,
    "rel_efficiency_cu",     "Relative efficiency",  "Complete",   2,
    "rel_efficiency_iu",     "Relative efficiency",  "Incomplete", 2,
    "bias_diff_complete",    "CATE bias difference", "All",        o$digits_bias,
    "bias_diff_complete_cu", "CATE bias difference", "Complete",   o$digits_bias,
    "bias_diff_complete_iu", "CATE bias difference", "Incomplete", o$digits_bias
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

caption_text <- function(o, model, runs) {
  paste0(
    model, ", missing-covariate study, ", o$dir, " outcome: mean (Monte Carlo ",
    "SE) over ", runs, " runs per cell. Both metrics are against the ",
    "complete-data arm, same run and estimator. Relative efficiency: MSE over ",
    "its MSE, 1 is no loss. CATE bias difference: CATE bias (mean of the ",
    "estimated minus true CATE) minus its CATE bias, 0 is no change. All: ",
    "every unit; complete: units with no amputed covariate; incomplete: the ",
    "rest. --- : not defined (complete case and IPW drop the incomplete units)."
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
  metrics <- readRDS(file.path(res_path, paste0(o$prefix, "_miss_metrics.RDS"))) %>%
    filter(scenario %in% SCENARIOS)
  runs <- runs_per_cell(metrics, run_cols)

  # --- estimation metrics -----------------------------------------------------
  cols <- table_cols(o)

  # the metrics file has the bias difference on all units and the complete
  # units only, so take the incomplete units' here, as thesis_figures/missing.R
  # does. Complete case and IPW's incomplete-unit bias is all NA, so theirs
  # comes out NA too.
  ref_bias_iu <- metrics %>%
    filter(method == REF_METHOD) %>%
    select(all_of(REF_BY), bias_iu_complete = bias_iu)

  metrics_summary <- metrics %>%
    filter(method != REF_METHOD) %>%
    left_join(ref_bias_iu, by = REF_BY) %>%
    mutate(bias_diff_complete_iu = bias_iu - bias_iu_complete) %>%
    summarise_metrics(run_cols, cols = setNames(cols$stem, cols$stem),
                      count_na = character()) %>%
    apply_labels(SS_SCENARIO_LABELS) %>%
    tex_mechanism()

  absent <- !paste0("mean_", cols$stem) %in% names(metrics_summary)
  if (any(absent)) {
    warning(o$prefix, "_miss_metrics.RDS has no ",
            paste(cols$stem[absent], collapse = ", "), " - columns dropped")
    cols <- cols[!absent, ]
  }

  # in MODEL_LABELS order; the file names keep the raw model name
  for (m in intersect(names(MODEL_LABELS), unique(metrics$model))) {
    tex <- metrics_summary %>%
      filter(model == MODEL_LABELS[[m]]) %>%
      droplevels() %>%
      miss_metrics_table(
        cols,
        caption.short = paste0(MODEL_LABELS[[m]], ", missing covariates, ",
                               o$dir, " outcome"),
        caption = caption_text(o, MODEL_LABELS[[m]], runs),
        label = paste0(o$prefix, "_miss_", m),
        # portrait fits the thesis's 17cm text width (A4, 2cm side margins) at
        # grouped_longtable()'s 8pt and 3pt padding: the widest binary row,
        # relative efficiency >= 10 in all three columns, is 16.8cm. At 9pt
        # the binary tables do not fit.
        landscape = FALSE
      )
    out_file <- file.path(tab_path, paste0(o$prefix, "_miss_", m, ".tex"))
    writeLines(tex, out_file)
    message("wrote ", out_file)
  }

  # --- HTE tests --------------------------------------------------------------
  true_cate_tests <- readRDS(
    file.path(res_path, paste0(o$prefix, "_miss_true_cate_tests.RDS"))
  ) %>%
    filter(scenario %in% SCENARIOS) %>%
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
    apply_labels(SS_SCENARIO_LABELS) %>%
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
