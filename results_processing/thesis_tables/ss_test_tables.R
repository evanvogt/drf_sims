##########
# title: LaTeX tables of the HTE tests - sample size, both outcomes
##########
# Rejection rates of the four heterogeneity tests, as "rate (MCSE)", with
# sample size across the columns so each test's power reads left to right.
# One table for the reported scenarios 1-4 and one for the supplementary 5-10,
# per outcome. The correlated-covariate studies (sample_size/correlated/) get
# the supplementary table only, once per rho: their scenarios 1-4 are the main
# text's figures (thesis_figures/sample_size.R). The estimation metrics are in
# ss_tables.R.
#
# Writes to ../results/thesis_tables/:
#   cts_ss_tests_main.tex, cts_ss_tests_supp.tex,
#   bin_ss_tests_main.tex, bin_ss_tests_supp.tex,
#   {cts,bin}_corr_rho0_ss_tests_supp.tex, {cts,bin}_corr_rho05_ss_tests_supp.tex
# \input{} them into a document with
#   \usepackage{booktabs, longtable, pdflscape, array}
#
# - Rejection is p < ALPHA, 0/1 per run, so the MCSE is the binomial one. Runs
#   with an NA p-value (a constant CATE estimate, or a test not run) are left
#   out of the denominator, as in the figures.
# - "True CATE" is each test run on the true CATE and nuisances
#   ({cts,bin}_true_cate_tests.RDS, from *_metrics.R): the ceiling the models'
#   own tests are chasing. It has no PO independence test (the DR-oracle's is
#   already the true pseudo-outcome's), and in scenario 1 the true CATE is
#   constant, so none of its tests are defined there.
# - A row with no defined cell is dropped: the PO independence test for the
#   causal forest and T-learners (they share another estimator's
#   pseudo-outcome, R/metrics.R), and the True CATE in scenario 1.

library(here)
source(here("R", "figures.R"))
source(here("R", "tables.R"))

ALPHA <- 0.05

# paths
path <- here()
tab_path <- file.path(dirname(path), "results", "thesis_tables")
dir.create(tab_path, showWarnings = FALSE, recursive = TRUE)

# Per entry: `dir` and `prefix` locate the metrics; `out` names the output
# files and LaTeX labels; `label` is the caption's outcome text; `rho` (the
# correlated studies only) is the rho kept; `tables` is which of main (1-4)
# and supp (5-10) to write; `note` is appended to the caption.
outcomes <- list(
  list(dir = "continuous", prefix = "cts", out = "cts",
       label = "continuous outcome", tables = c("main", "supp"), note = ""),
  list(dir = "binary", prefix = "bin", out = "bin",
       label = "binary outcome", tables = c("main", "supp"), note = "")
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

# p-value column -> row label, in row order
tests <- c(
  BLP_p = "BLP (two-sided)",
  BLP_p_os = "BLP (one-sided, HC3)",
  indep_cate = "CATE independence",
  indep_po = "PO independence"
)
rej_tests <- setNames(tests, paste0("rej_", names(tests)))

# p-values -> 0/1 rejection indicators, rej_<test>. A test absent from `df`
# (indep_po in the true-CATE file) is left absent.
add_rejections <- function(df) {
  df %>%
    mutate(across(any_of(names(tests)), ~ as.numeric(as.numeric(.x) < ALPHA),
                  .names = "rej_{.col}"))
}

summarise_rejections <- function(df, group_cols) {
  summarise_metrics(
    df,
    group_cols,
    cols = setNames(names(rej_tests), names(rej_tests)),
    binomial = names(rej_tests),
    count_na = character()
  )
}

short_caption_text <- function(o, scenarios) {
  paste0("HTE tests, ", o$label, ", scenarios ", scenarios)
}

caption_text <- function(o, scenarios, runs) {
  paste0(
    "HTE tests, sample-size study, ", o$label, ", scenarios ", scenarios,
    ": rejection rate at $\\alpha = ", ALPHA, "$ (Monte Carlo SE) over ", runs,
    " runs per cell, by sample size. Type I error in the null scenario, power ",
    "otherwise. Runs with an undefined p-value are left out. True CATE: the ",
    "test applied to the true CATE and nuisances.", o$note
  )
}

# the rho this entry keeps, if it is a correlated study
keep_rho <- function(df, o) if (is.null(o$rho)) df else filter(df, rho == o$rho)

for (o in outcomes) {
  res_path <- file.path(dirname(path), "results", o$dir)
  metrics <- readRDS(file.path(res_path, paste0(o$prefix, "_metrics.RDS"))) %>%
    keep_rho(o) %>%
    add_rejections()
  true_cate_tests <- readRDS(
    file.path(res_path, paste0(o$prefix, "_true_cate_tests.RDS"))
  ) %>%
    keep_rho(o) %>%
    add_rejections()

  test_summary <- bind_rows(
    summarise_rejections(metrics, c("scenario", "n", "model")),
    summarise_rejections(true_cate_tests, c("scenario", "n")) %>%
      mutate(model = "True CATE")
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
    # apply_labels() puts the unlabelled "True CATE" after the estimators
    tex <- test_summary %>%
      filter(scenario %in% t$scenarios) %>%
      apply_labels(t$labels) %>%
      mutate(n = factor(n, levels = c(100, 250, 500, 1000))) %>%
      ss_test_table(
        rej_tests,
        digits = 3,
        caption.short = short_caption_text(o, t$text),
        caption = caption_text(
          o, t$text, runs_per_cell(filter(metrics, scenario %in% t$scenarios))
        ),
        label = paste0(o$out, "_ss_tests_", nm),
        landscape = FALSE
      )
    out_file <- file.path(tab_path, paste0(o$out, "_ss_tests_", nm, ".tex"))
    writeLines(tex, out_file)
    message("wrote ", out_file)
  }
}
