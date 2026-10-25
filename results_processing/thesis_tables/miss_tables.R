##########
# title: LaTeX tables for the thesis appendix - missing covariates, both outcomes
##########
# Appendix tables for the missing-data studies, scenarios 1-4, for the
# estimators in KEEP_MODELS: one estimation table per estimator (relative
# efficiency and CATE bias difference against the complete-data arm, same run,
# on the complete and the incomplete units, the three mechanisms across the
# columns) and one per HTE test. Cells are the mean (estimation) or the
# rejection rate (tests) alone, with the Monte Carlo SEs in the caption.
# Measures as in missing/ADEMP.md, "Performance measures".
#
# - The complete-data arm is split by the amputation its run would have had, so
#   each unit set is compared with the same units.
# - Complete case and IPW drop the incomplete units, so they have complete-unit
#   columns only.
# - Bias is the CATE bias (mean of est - true over units).
# - No "all units" columns: they mix the complete and incomplete ones, which
#   the chapter reads separately.
# - No relative bias (it blows up where the true CATE is near 0) and no Pearson
#   (the complete units' narrower CATE range lowers it even with no missingness).
# - The truth is tau(X), under MNAR-tau too.
# - Only the one-sided BLP and PO-independence tests are tabulated, as in the
#   sample-size appendix.
#
# Writes to ../results/thesis_tables/, for prefix cts and bin:
# <prefix>_miss_<model>.tex and <prefix>_miss_tests_<test>.tex. \input{} with
# \usepackage{booktabs, longtable, pdflscape, array}.
#
#   cd results_processing/thesis_tables; Rscript miss_tables.R
#   Rscript miss_tables.R [results root] [table dir]   # mock smoke run

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

# the estimators tabulated. The semi-oracle behaves as the DR-RF, and the
# oracle and T-learners are references in the main text, not appendix tables
KEEP_MODELS <- c("dr_random_forest", "causal_forest", "dr_superlearner")

# paths
path <- here()
args <- commandArgs(trailingOnly = TRUE)
res_root <- if (length(args) >= 1) args[1] else file.path(dirname(path), "results")
tab_path <- if (length(args) >= 2) args[2] else file.path(res_root, "thesis_tables")
dir.create(tab_path, showWarnings = FALSE, recursive = TRUE)

# the binary outcome's errors are on the risk-difference scale, an order of
# magnitude smaller, so bias gets an extra decimal
outcomes <- list(
  list(dir = "continuous", prefix = "cts", digits_bias = 3),
  list(dir = "binary", prefix = "bin", digits_bias = 4)
)

# the estimation tables' columns, left to right, under a spanning header per
# mechanism and per metric. Drop a row here to drop a column.
table_cols <- function(o) {
  tibble::tribble(
    ~stem,                   ~metric,                ~header,      ~digits,
    "rel_efficiency_cu",     "Rel. eff.",             "CU", 2,
    "rel_efficiency_iu",     "Rel. eff.",             "IU", 2,
    "bias_diff_complete_cu", "$\\\\Delta$ CATE bias", "CU", o$digits_bias,
    "bias_diff_complete_iu", "$\\\\Delta$ CATE bias", "IU", o$digits_bias
  )
}

# p-value column -> table title, in output order
tests <- c(
  BLP_p_os = "BLP (one-sided, HC3)",
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
# typeset. The mechanism heads a spanning header, where add_header_above()
# eats one backslash, so it takes two.
tex_mechanism <- function(df) {
  lv <- levels(df$mechanism)
  levels(df$mechanism)[lv == MECHANISM_LABELS[["MNAR-tau"]]] <- "MNAR-$\\\\tau$"
  df
}

group_cols <- c("scenario", "mechanism", "method")
run_cols <- c(group_cols, "model")

caption_text <- function(o, model, runs, mcse) {
  paste0(
    model, ", missing-covariate study, ", o$dir, " outcome: mean over ", runs,
    " runs per cell, by missingness mechanism. ", mcse, " Both metrics are ",
    "against the complete-data arm, same run and estimator. Rel. eff.: ",
    "relative efficiency, MSE over its MSE, 1 is no loss. $\\Delta$ CATE ",
    "bias: CATE bias difference, CATE bias (mean of the estimated minus true ",
    "CATE) minus its CATE bias, 0 is no change. CU: complete units, with no ",
    "amputed covariate; IU: incomplete units, the rest. MIM: missing ",
    "indicator method. --- : not defined (complete case ",
    "and IPW drop the incomplete units)."
  )
}

test_caption_text <- function(o, test, runs, mcse) {
  paste0(
    test, " test, missing-covariate study, ", o$dir, " outcome: rejection ",
    "rate at $\\alpha = ", ALPHA, "$ over ", runs, " runs per cell, by ",
    "missingness mechanism. ", mcse, " Type I error in the null scenario, ",
    "power otherwise. Runs with an undefined p-value are left out, and ",
    "handling methods with none at all. True CATE: the test applied to the ",
    "true CATE and nuisances, at the handled covariates. MIM: missing ",
    "indicator method. --- : not run or not defined."
  )
}

for (o in outcomes) {
  res_path <- file.path(res_root, "missing", o$dir)
  metrics <- readRDS(file.path(res_path, paste0(o$prefix, "_miss_metrics.RDS"))) %>%
    filter(scenario %in% SCENARIOS)
  runs <- runs_per_cell(metrics, run_cols)
  run_counts <- count(metrics, across(all_of(run_cols)), name = "runs")$runs

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
  for (m in intersect(intersect(names(MODEL_LABELS), KEEP_MODELS),
                      unique(metrics$model))) {
    summ <- metrics_summary %>%
      filter(model == MODEL_LABELS[[m]]) %>%
      droplevels()

    d <- bind_rows(lapply(seq_len(nrow(cols)), function(i) {
      transmute(
        summ,
        group = scenario, method, mechanism,
        metric = factor(cols$metric[i], levels = unique(cols$metric)),
        col = factor(cols$header[i], levels = unique(cols$header)),
        cell = fmt_num(.data[[paste0("mean_", cols$stem[i])]], cols$digits[i])
      )
    }))

    # one SE per metric, over its complete and incomplete columns
    worst <- summ %>%
      mutate(
        mcse_re = pmax(mcse_rel_efficiency_cu, mcse_rel_efficiency_iu, na.rm = TRUE),
        mcse_bd = pmax(mcse_bias_diff_complete_cu, mcse_bias_diff_complete_iu,
                       na.rm = TRUE)
      )

    tex <- wide_table(
      d,
      row_cols = c(Method = "method"),
      spans = c("metric", "mechanism"),
      col = "col",
      caption.short = paste0(MODEL_LABELS[[m]], ", missing covariates, ",
                             o$dir, " outcome"),
      caption = caption_text(
        o, MODEL_LABELS[[m]], runs,
        mcse_note(worst, c(re = "relative efficiency", bd = "CATE bias difference"),
                  c(2, o$digits_bias))
      ),
      label = paste0(o$prefix, "_miss_", m),
      # 12 numeric columns at 8pt and 3pt padding fit the thesis's 17cm text
      # width (A4, 2cm side margins)
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
    summarise_rejections(
      add_rejections(filter(metrics, model %in% KEEP_MODELS)), run_cols
    ),
    summarise_rejections(true_cate_tests, group_cols) %>%
      mutate(model = "True CATE")
  ) %>%
    apply_labels(SS_SCENARIO_LABELS) %>%
    tex_mechanism()

  for (test in names(tests)) {
    stem <- paste0("rej_", test)
    cell <- fmt_num(test_summary[[paste0("mean_", stem)]], 3)
    d <- test_summary %>%
      transmute(group = scenario, method, mechanism, model, cell = cell) %>%
      # a column with no defined cell (the causal forest's PO independence)
      group_by(model) %>%
      filter(any(cell != "---")) %>%
      ungroup()

    tex <- wide_table(
      d,
      row_cols = c(Method = "method"),
      spans = "mechanism",
      col = "model",
      caption.short = paste0(tests[[test]], " test, missing covariates, ",
                             o$dir, " outcome"),
      caption = test_caption_text(o, tests[[test]], runs,
                                  binom_mcse_note(unique(range(run_counts)))),
      label = paste0(o$prefix, "_miss_tests_", test),
      landscape = FALSE
    )
    out_file <- file.path(tab_path, paste0(o$prefix, "_miss_tests_", test, ".tex"))
    writeLines(tex, out_file)
    message("wrote ", out_file)
  }
}
