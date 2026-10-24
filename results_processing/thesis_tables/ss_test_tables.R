##########
# title: LaTeX tables of the HTE tests - sample size, both outcomes
##########
# Rejection rates of the one-sided BLP test and the PO-independence test, with
# n across the columns, for the correlated sample-size studies: rho = 0.5 and
# the paired difference rho = 0.5 minus rho = 0, scenarios 1-4 (main) and 5-10
# (supp). Continuous and binary outcome sit side by side in one table, cells
# the rate alone, with the Monte Carlo SE in the caption. The two-sided BLP
# test and the CATE-independence test are summarised in the chapter's text and
# figures, not tabulated.
#
# - Rejection is p < HTE_ALPHA, 0/1 per run (binomial MCSE sqrt(p(1-p)/R)). NA
#   p-values are left out of the denominator, and out of a pair.
# - "True CATE" is each test on the true CATE and nuisances: the ceiling the
#   models chase. It has no PO independence test and is undefined in scenario 1.
# - Rows with no defined cell are dropped.
#
# Writes to ../results/thesis_tables/:
#   corr_{rho05,rhodiff}_ss_tests_{main,supp}.tex
# \input{} with \usepackage{booktabs, longtable, pdflscape, array}.
#
#   cd results_processing/thesis_tables; Rscript ss_test_tables.R
#   Rscript ss_test_tables.R [results root] [table dir]   # mock smoke run

library(here)
source(here("R", "figures.R"))
source(here("R", "tables.R"))

# paths
path <- here()
args <- commandArgs(trailingOnly = TRUE)
res_root <- if (length(args) >= 1) args[1] else file.path(dirname(path), "results")
tab_path <- if (length(args) >= 2) args[2] else file.path(res_root, "thesis_tables")
dir.create(tab_path, showWarnings = FALSE, recursive = TRUE)

# `supp_note` is appended to the 5-10 table's caption: binary scenario 9 runs
# on a smaller scale than binary/'s (RD_SCALE_CORR in R/dgm_scenarios.R)
outcomes <- list(
  list(dir = "continuous", prefix = "cts", header = "Continuous"),
  list(dir = "binary", prefix = "bin", header = "Binary")
)
SUPP_NOTE <- paste0(
  " In binary scenario 9 the heterogeneity is scaled by 0.139 rather than the",
  " independent-covariate study's 0.164, so that every treated risk",
  " stays within $[0.01, 0.99]$ at $\\rho = 0.5$."
)

# the table sets: the rho = 0.5 rates and the paired differences (CORR_RHOS in
# R/dgm_scenarios.R). `out` names the output files and LaTeX labels, `label`
# is the caption's covariate text
sets <- list(
  list(out = "corr_rho05", label = "correlated covariates ($\\rho = 0.5$)",
       rho = 0.5, diff = FALSE),
  list(out = "corr_rhodiff",
       label = paste0("correlated covariates ($\\rho = 0.5$) minus independent ",
                      "covariates ($\\rho = 0$)"),
       diff = TRUE)
)

# p-value column -> row label, in row order
tests <- c(
  BLP_p_os = "BLP (one-sided, HC3)",
  indep_po = "PO independence"
)
rej_tests <- setNames(tests, paste0("rej_", names(tests)))

# p-values -> 0/1 rejection indicators, rej_<test>. A test absent from `df`
# (indep_po in the true-CATE file) is left absent.
add_rejections <- function(df) {
  df %>%
    mutate(across(any_of(names(tests)),
                  ~ as.numeric(as.numeric(.x) < HTE_ALPHA),
                  .names = "rej_{.col}"))
}

#' Rates get the binomial MCSE; differences (diff = TRUE) the general one,
#' which on per-run differences is the paired SE
summarise_rejections <- function(df, group_cols, diff) {
  summarise_metrics(
    df,
    group_cols,
    cols = setNames(names(rej_tests), names(rej_tests)),
    binomial = if (diff) character() else names(rej_tests),
    count_na = character()
  )
}

short_caption_text <- function(s, scenarios) {
  paste0("HTE tests, ", s$label, ", scenarios ", scenarios)
}

caption_text <- function(s, scenarios, runs, mcse, note = "") {
  cells <- if (s$diff) {
    paste0(
      "per-run difference in rejection at $\\alpha = ", HTE_ALPHA, "$, by ",
      "sample size, over ", runs, " run pairs per cell. A run shares its ",
      "random draws at the two values of $\\rho$, so the Monte Carlo SE is ",
      "that of the paired difference. "
    )
  } else {
    paste0(
      "rejection rate at $\\alpha = ", HTE_ALPHA, "$, by sample size, over ",
      runs, " runs per cell. Type I error in the null scenario, power ",
      "otherwise. "
    )
  }
  paste0(
    "HTE tests, sample-size study, ", s$label, ", scenarios ", scenarios,
    ": ", cells, mcse, " Runs with an undefined p-value are left out. ",
    "True CATE: the test applied to the true CATE and nuisances.", note
  )
}

for (s in sets) {
  # both outcomes' rows, one `outcome` column apart
  per_outcome <- lapply(outcomes, function(o) {
    res_path <- file.path(res_root, "correlated", o$dir)
    # both rhos, then the rows at s$rho or the per-run differences
    read_corr <- function(stem, keys) {
      df <- readRDS(file.path(res_path, paste0(o$prefix, "_corr_", stem, ".RDS"))) %>%
        add_rejections()
      if (s$diff) {
        paired_rho_diff(df, intersect(names(rej_tests), names(df)), keys)
      } else {
        filter(df, rho == s$rho)
      }
    }
    metrics <- read_corr("metrics", c("scenario", "n", "model", "run"))
    true_cate_tests <- read_corr("true_cate_tests", c("scenario", "n", "run"))

    list(
      metrics = metrics,
      summary = bind_rows(
        summarise_rejections(metrics, c("scenario", "n", "model"), s$diff),
        summarise_rejections(true_cate_tests, c("scenario", "n"), s$diff) %>%
          mutate(model = "True CATE")
      ) %>%
        mutate(outcome = o$header)
    )
  })
  names(per_outcome) <- vapply(outcomes, `[[`, "", "header")
  metrics <- bind_rows(lapply(per_outcome, `[[`, "metrics"), .id = "outcome")
  test_summary <- bind_rows(lapply(per_outcome, `[[`, "summary")) %>%
    mutate(outcome = factor(outcome, levels = vapply(outcomes, `[[`, "", "header")))

  tables <- list(
    main = list(scenarios = 1:4, labels = SS_SCENARIO_LABELS,
                text = "1--4", note = ""),
    supp = list(scenarios = 5:10, labels = SS_SUPP_SCENARIO_LABELS,
                text = "5--10", note = SUPP_NOTE)
  )

  for (nm in names(tables)) {
    t <- tables[[nm]]
    # 5-10 were added to the studies after 1-4 had run
    if (!any(metrics$scenario %in% t$scenarios)) {
      message(s$out, ": no scenario ", t$text, " metrics yet - ", nm, " table skipped")
      next
    }
    # apply_labels() puts the unlabelled "True CATE" after the estimators
    summ <- test_summary %>%
      filter(scenario %in% t$scenarios) %>%
      apply_labels(t$labels) %>%
      mutate(n = factor(n, levels = c(100, 250, 500, 1000)))

    d <- bind_rows(lapply(names(rej_tests), function(stem) {
      transmute(
        summ,
        group = scenario, test = factor(rej_tests[[stem]], levels = unname(tests)),
        model, outcome, n,
        cell = fmt_num(.data[[paste0("mean_", stem)]], 3)
      )
    }))

    in_set <- filter(metrics, scenario %in% t$scenarios)
    runs <- count(in_set, scenario, n, model, outcome, name = "runs")$runs
    mcse <- if (s$diff) {
      summ %>%
        mutate(mcse_diff = pmax(!!!syms(paste0("mcse_", names(rej_tests))),
                                na.rm = TRUE)) %>%
        mcse_note(c(diff = "difference in rejection rate"), 3)
    } else {
      binom_mcse_note(unique(range(runs)))
    }

    tex <- wide_table(
      d,
      row_cols = c(Test = "test", Model = "model"),
      spans = "outcome",
      col = "n",
      block = TRUE,
      caption.short = short_caption_text(s, t$text),
      caption = caption_text(
        s, t$text, runs_per_cell(in_set, c("scenario", "n", "model", "outcome")),
        mcse, t$note
      ),
      label = paste0(s$out, "_ss_tests_", nm),
      landscape = FALSE
    )
    out_file <- file.path(tab_path, paste0(s$out, "_ss_tests_", nm, ".tex"))
    writeLines(tex, out_file)
    message("wrote ", out_file)
  }
}
