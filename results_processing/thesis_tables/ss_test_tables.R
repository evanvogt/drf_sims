##########
# title: LaTeX tables of the HTE tests - sample size, both outcomes
##########
# Rejection rates of the four heterogeneity tests, as "rate (MCSE)", with
# sample size across the columns so each test's power reads left to right.
# One table for the reported scenarios 1-4 and one for the supplementary 5-10,
# per outcome. The estimation metrics are in ss_tables.R.
#
# Writes to ../results/thesis_tables/:
#   cts_ss_tests_main.tex, cts_ss_tests_supp.tex,
#   bin_ss_tests_main.tex, bin_ss_tests_supp.tex
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

outcomes <- list(
  list(dir = "continuous", prefix = "cts"),
  list(dir = "binary", prefix = "bin")
)

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
  paste0("HTE tests, ", o$dir, " outcome, scenarios ", scenarios)
}

caption_text <- function(o, scenarios, runs) {
  paste0(
    "HTE tests, sample-size study, ", o$dir, " outcome, scenarios ", scenarios,
    ": rejection rate at $\\alpha = ", ALPHA, "$ (Monte Carlo SE) over ", runs,
    " runs per cell, by sample size. Type I error in the null scenario, power ",
    "otherwise. Runs with an undefined p-value are left out. True CATE: the ",
    "test applied to the true CATE and nuisances."
  )
}

for (o in outcomes) {
  res_path <- file.path(dirname(path), "results", o$dir)
  metrics <- readRDS(file.path(res_path, paste0(o$prefix, "_metrics.RDS"))) %>%
    add_rejections()
  true_cate_tests <- readRDS(
    file.path(res_path, paste0(o$prefix, "_true_cate_tests.RDS"))
  ) %>%
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

  for (nm in names(tables)) {
    t <- tables[[nm]]
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
        label = paste0(o$prefix, "_ss_tests_", nm),
        landscape = FALSE
      )
    out_file <- file.path(tab_path, paste0(o$prefix, "_ss_tests_", nm, ".tex"))
    writeLines(tex, out_file)
    message("wrote ", out_file)
  }
}
