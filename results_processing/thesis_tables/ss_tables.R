##########
# title: LaTeX tables for the thesis chapter - sample size, both outcomes
##########
# Appendix tables for the correlated sample-size studies, per outcome: rho = 0.5
# and the paired difference rho = 0.5 minus rho = 0 (the rho = 0 level is the
# first minus the second, so it has no table of its own), each for scenarios
# 1-4 (main) and 5-10 (supp). Layout: one row per scenario x model, the four
# sample sizes across the columns under each metric, cells the mean alone. The
# largest Monte Carlo SE of each metric is in the caption. The HTE tests are in
# ss_test_tables.R.
#
# - Runs share draws across rho, so differences are per run, MCSE
#   sd(diff) / sqrt(pairs); Pearson gets a third decimal there.
# - Bias is the CATE bias (mean of est - true), not an ATE estimate.
# - No relative bias: the true ATE is recalibrated per n, so it has no common
#   scale across n.
# - Pearson, not Spearman: the true CATE is mostly ties in scenarios 5 and 7.
#   Pearson is undefined in scenario 1 and prints as a dash.
# - No sign accuracy: the chapter's results never report it.
#
# Writes to ../results/thesis_tables/:
#   {cts,bin}_corr_{rho05,rhodiff}_ss_{main,supp}.tex
# \input{} with \usepackage{booktabs, longtable, pdflscape, array}.
#
#   cd results_processing/thesis_tables; Rscript ss_tables.R
#   Rscript ss_tables.R [results root] [table dir]   # mock smoke run

library(here)
source(here("R", "figures.R"))
source(here("R", "tables.R"))

# paths
path <- here()
args <- commandArgs(trailingOnly = TRUE)
res_root <- if (length(args) >= 1) args[1] else file.path(dirname(path), "results")
tab_path <- if (length(args) >= 2) args[2] else file.path(res_root, "thesis_tables")
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

# one entry per outcome and table set, the rho = 0.5 level and the paired
# difference (CORR_RHOS in R/dgm_scenarios.R): `out` names the output files and
# LaTeX labels, `label` is the caption's outcome text
outcomes <- unlist(lapply(outcomes, function(o) list(
  modifyList(o, list(
    out = paste0(o$prefix, "_corr_rho05"),
    label = paste0(o$label, ", correlated covariates ($\\rho = 0.5$)"),
    rho = 0.5,
    diff = FALSE
  )),
  modifyList(o, list(
    out = paste0(o$prefix, "_corr_rhodiff"),
    label = paste0(o$label, ", correlated covariates ($\\rho = 0.5$) minus ",
                   "independent covariates ($\\rho = 0$)"),
    diff = TRUE
  ))
)), recursive = FALSE)

# the table's metrics, left to right. Drop a row here to drop a metric.
table_cols <- function(o) {
  digits_prop <- if (o$diff) 3 else 2
  tibble::tribble(
    ~stem,  ~header,    ~name,      ~digits,
    "bias", "Bias",     "bias",     o$digits_bias,
    "mse",  "MSE",      "MSE",      o$digits_mse,
    "corr", "Pearson",  "Pearson correlation", digits_prop
  )
}

short_caption_text <- function(o, scenarios) {
  paste0("CATE estimation, ", o$label, ", scenarios ", scenarios)
}

caption_text <- function(o, scenarios, runs, mcse, note = "") {
  cells <- if (o$diff) {
    paste0(
      "per-run difference, mean over ", runs, " run pairs per cell. A run ",
      "shares its random draws at the two values of $\\rho$, so the Monte ",
      "Carlo SE is that of the paired difference. "
    )
  } else {
    paste0("mean over ", runs, " runs per cell. ")
  }
  paste0(
    "CATE estimation, sample-size study, ", o$label, ", scenarios ",
    scenarios, ": ", cells, mcse, " Bias is the CATE bias, the mean of the ",
    "estimated minus true CATE over units. --- : not defined (correlation ",
    "under no heterogeneity).", note
  )
}

# paired_rho_diff() is R/figures.R's, shared with ss_test_tables.R and the
# figures
for (o in outcomes) {
  metrics <- readRDS(file.path(res_root, "correlated", o$dir,
                               paste0(o$prefix, "_corr_metrics.RDS")))

  metrics <- metrics %>%
    mutate(corr = if_else(scenario == 1, NA_real_, corr))

  cols <- table_cols(o)

  metrics <- if (o$diff) {
    paired_rho_diff(metrics, cols$stem)
  } else {
    filter(metrics, rho == o$rho)
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
    summ <- metrics_summary %>%
      filter(scenario %in% t$scenarios) %>%
      apply_labels(t$labels) %>%
      mutate(n = factor(n, levels = c(100, 250, 500, 1000)))

    # the cells where the SE is not small against the mean: the caption's
    # maximum is only a fair summary if these are few
    for (i in which(cols$stem != "bias" & !o$diff)) {
      s <- cols$stem[i]
      noisy <- summ[[paste0("mcse_", s)]] > 0.3 * abs(summ[[paste0("mean_", s)]])
      if (any(noisy, na.rm = TRUE)) {
        message(o$out, "_ss_", nm, ": ", sum(noisy, na.rm = TRUE), " of ",
                nrow(summ), " ", cols$header[i], " cells have an MCSE over 30% ",
                "of the mean:")
        print(as.data.frame(summ[which(noisy), c("scenario", "n", "model")]) %>%
                mutate(mean = summ[[paste0("mean_", s)]][which(noisy)],
                       mcse = summ[[paste0("mcse_", s)]][which(noisy)]),
              digits = 3)
      }
    }

    d <- bind_rows(lapply(seq_len(nrow(cols)), function(i) {
      transmute(
        summ,
        group = scenario, model, n,
        metric = factor(cols$header[i], levels = cols$header),
        cell = fmt_num(.data[[paste0("mean_", cols$stem[i])]], cols$digits[i])
      )
    }))

    tex <- wide_table(
      d,
      row_cols = c(Model = "model"),
      spans = "metric",
      col = "n",
      caption.short = short_caption_text(o, t$text),
      caption = caption_text(
        o, t$text, runs_per_cell(filter(metrics, scenario %in% t$scenarios)),
        mcse_note(summ, setNames(cols$name, cols$stem), cols$digits),
        t$note
      ),
      label = paste0(o$out, "_ss_", nm),
      # 12 numeric columns: 8pt, as every other appendix table
      landscape = FALSE
    )
    out_file <- file.path(tab_path, paste0(o$out, "_ss_", nm, ".tex"))
    writeLines(tex, out_file)
    message("wrote ", out_file)
  }
}
