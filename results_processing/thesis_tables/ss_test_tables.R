##########
# title: LaTeX tables of the HTE tests - sample size, both outcomes
##########
# Rejection rates of the four heterogeneity tests, as "rate (MCSE)", with
# sample size across the columns so each test's power reads left to right.
# The correlated-covariate studies (sample_size/correlated/), per outcome and
# rho: one table for the reported scenarios 1-4 and one for the supplementary
# 5-10. The independent-covariate studies (sample_size/continuous/, binary/)
# are no longer tabulated. The estimation metrics are in ss_tables.R.
#
# Writes to ../results/thesis_tables/:
#   {cts,bin}_corr_rho0_ss_tests_main.tex, {cts,bin}_corr_rho0_ss_tests_supp.tex,
#   {cts,bin}_corr_rho05_ss_tests_main.tex, {cts,bin}_corr_rho05_ss_tests_supp.tex
# \input{} them into a document with
#   \usepackage{booktabs, longtable, pdflscape, array}
#
# - Rejection is p < ALPHA, 0/1 per run, so the MCSE is the binomial one. Runs
#   with an NA p-value (a constant CATE estimate, or a test not run) are left
#   out of the denominator, as in the figures.
# - "True CATE" is each test run on the true CATE and nuisances
#   ({cts,bin}_corr_true_cate_tests.RDS, from *_corr_metrics.R): the ceiling
#   the models' own tests are chasing. It has no PO independence test (the
#   DR-oracle's is already the true pseudo-outcome's), and in scenario 1 the
#   true CATE is constant, so none of its tests are defined there.
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

# `supp_note` is appended to the 5-10 table's caption: binary scenario 9 runs
# on a smaller scale than binary/'s (RD_SCALE_CORR in R/dgm_scenarios.R)
outcomes <- list(
  list(dir = "continuous", prefix = "cts", label = "continuous outcome",
       supp_note = ""),
  list(dir = "binary", prefix = "bin", label = "binary outcome",
       supp_note = paste0(
         " In scenario 9 the heterogeneity is scaled by 0.139 rather than the",
         " independent-covariate study's 0.164, so that every treated risk",
         " stays within $[0.01, 0.99]$ at $\\rho = 0.5$."
       ))
)

# one entry per outcome and rho (CORR_RHOS in R/dgm_scenarios.R): `out` names
# the output files and LaTeX labels, `label` is the caption's outcome text
outcomes <- unlist(lapply(outcomes, function(o) {
  lapply(c(0, 0.5), function(r) modifyList(o, list(
    out = paste0(o$prefix, "_corr_rho", sub(".", "", r, fixed = TRUE)),
    label = paste0(o$label, ", correlated covariates ($\\rho = ", r, "$)"),
    rho = r
  )))
}), recursive = FALSE)

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

caption_text <- function(o, scenarios, runs, note = "") {
  paste0(
    "HTE tests, sample-size study, ", o$label, ", scenarios ", scenarios,
    ": rejection rate at $\\alpha = ", ALPHA, "$ (Monte Carlo SE) over ", runs,
    " runs per cell, by sample size. Type I error in the null scenario, power ",
    "otherwise. Runs with an undefined p-value are left out. True CATE: the ",
    "test applied to the true CATE and nuisances.", note
  )
}

for (o in outcomes) {
  res_path <- file.path(dirname(path), "results", "correlated", o$dir)
  read_corr <- function(stem) {
    readRDS(file.path(res_path, paste0(o$prefix, "_corr_", stem, ".RDS"))) %>%
      filter(rho == o$rho) %>%
      add_rejections()
  }
  metrics <- read_corr("metrics")
  true_cate_tests <- read_corr("true_cate_tests")

  test_summary <- bind_rows(
    summarise_rejections(metrics, c("scenario", "n", "model")),
    summarise_rejections(true_cate_tests, c("scenario", "n")) %>%
      mutate(model = "True CATE")
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
          o, t$text, runs_per_cell(filter(metrics, scenario %in% t$scenarios)),
          t$note
        ),
        label = paste0(o$out, "_ss_tests_", nm),
        landscape = FALSE
      )
    out_file <- file.path(tab_path, paste0(o$out, "_ss_tests_", nm, ".tex"))
    writeLines(tex, out_file)
    message("wrote ", out_file)
  }
}
