##########
# title: LaTeX tables of the CI studies - sample size, both outcomes
##########
# Appendix tables for the correlated CI studies (sample_size/confidence_intervals/):
# the CI_sf sweep per outcome (rho = 0, rho = 0.5 and their paired difference)
# and the data-driven CI_sf (optimal_sf/, continuous), plus the main-text float
# tab:optimal_sf. Cells are "mean (MCSE)".
#
# - Per-unit intervals only: the query grid is extrapolation at rho = 0.5
#   (sample_size/ADEMP.md, "Estimands").
# - causal_forest_inbuilt never sees CI_sf, so it is printed once across the
#   ratio columns.
# - Simultaneous coverage is 0/1 per run (binomial MCSE); marginal coverage is a
#   within-run proportion (general MCSE). Differences are per run, MCSE
#   sd(diff) / sqrt(pairs).
# - optimal_sf bands are 90%, the sweep's 95%.
# - A study with missing or incomplete metrics is skipped with a message.
#
# Writes to ../results/thesis_tables/:
#   {cts,bin}_corr_{rho0,rho05,rhodiff}_ci_{marg,simul,len}.tex
#   cts_corr_ci_sf.tex, cts_corr_ci_sf_main.tex
# \input{} with \usepackage{booktabs, longtable, pdflscape, array}.

library(here)
source(here("R", "figures.R"))
source(here("R", "tables.R"))

# paths
path <- here()
res_root <- file.path(dirname(path), "results", "correlated", "confidence_intervals")
tab_path <- file.path(dirname(path), "results", "thesis_tables")
dir.create(tab_path, showWarnings = FALSE, recursive = TRUE)

# the studies' designs (*_corr_ci_config.R, *_corr_ci_sf_config.R), to check
# the metrics are complete before tabulating them
DESIGN <- expand.grid(rho = c(0, 0.5), scenario = 1:4, n = c(500, 1000))
CI_SF <- round(seq(0.05, 0.5, 0.05), 2)

# the binary outcome's intervals are on the risk-difference scale, an order of
# magnitude narrower, so length gets an extra decimal
outcomes <- list(
  list(dir = "continuous", prefix = "cts", label = "continuous outcome",
       digits_len = 2),
  list(dir = "binary", prefix = "bin", label = "binary outcome",
       digits_len = 3)
)

# the sweep's metrics, one table each: `file` names the output, `binomial`
# marks a per-run 0/1 indicator. `desc` opens the caption, `defn` explains
# the metric in it.
ci_metrics <- tibble::tribble(
  ~stem,       ~col,                    ~file,   ~digits, ~binomial,
  "marg_cov",  "marginal_coverage",     "marg",  3,       FALSE,
  "simul_cov", "simultaneous_coverage", "simul", 2,       TRUE,
  "ci_len",    "mean_ci_length",        "len",   NA,      FALSE
) %>%
  mutate(
    desc = c("Marginal coverage", "Simultaneous coverage", "Mean interval length"),
    defn = c(
      paste0("Marginal coverage is the share of units whose true CATE the ",
             "interval covers. Nominal coverage is 0.95. "),
      paste0("Simultaneous coverage is the share of runs in which the interval ",
             "covers every unit's true CATE at once. Nominal coverage is 0.95. "),
      "Length is averaged over units. "
    )
  )

# the one interval that does not depend on CI_sf (see the header)
SPAN_MODELS <- MODEL_LABELS[["causal_forest_inbuilt"]]

# one entry per outcome and rho table set, as ss_tables.R: `out` names the
# output files and LaTeX labels, `label` is the caption's outcome text
sets <- unlist(lapply(outcomes, function(o) list(
  modifyList(o, list(
    out = paste0(o$prefix, "_corr_rho0"),
    label = paste0(o$label, ", correlated covariates ($\\rho = 0$)"),
    rho = 0,
    diff = FALSE
  )),
  modifyList(o, list(
    out = paste0(o$prefix, "_corr_rho05"),
    label = paste0(o$label, ", correlated covariates ($\\rho = 0.5$)"),
    rho = 0.5,
    diff = FALSE
  )),
  modifyList(o, list(
    out = paste0(o$prefix, "_corr_rhodiff"),
    label = paste0(o$label, ", correlated covariates, $\\rho = 0.5$ minus ",
                   "$\\rho = 0$"),
    diff = TRUE
  ))
)), recursive = FALSE)

#' The design cells (rows of `design`, columns `by`) with no row in `metrics`
missing_cells <- function(metrics, design, by) {
  anti_join(design, distinct(metrics, across(all_of(by))), by = by)
}

#' Read a metrics file, or NULL with a message if it is missing or incomplete
read_complete <- function(file, design, by, what) {
  if (!file.exists(file)) {
    message(what, ": ", file, " not found - skipped")
    return(NULL)
  }
  metrics <- readRDS(file)
  # CI_sf is stored as seq()'s doubles (0.15 is 0.15000000000000002), so it
  # is rounded here and in CI_SF, and the join is exact
  if ("CI_sf" %in% names(metrics)) metrics <- mutate(metrics, CI_sf = round(CI_sf, 2))
  gap <- missing_cells(metrics, design, by)
  if (nrow(gap) > 0) {
    message(what, ": ", nrow(gap), " of ", nrow(design), " design cells have ",
            "no metrics yet - skipped")
    return(NULL)
  }
  metrics
}

sweep_caption <- function(o, m, runs) {
  cells <- if (o$diff) {
    paste0(
      "per-run difference, mean (Monte Carlo SE) over ", runs, " run pairs ",
      "per cell. A run shares its random draws at the two values of $\\rho$, ",
      "so the SE is that of the paired difference. "
    )
  } else {
    paste0("mean (Monte Carlo SE) over ", runs, " runs per cell. ")
  }
  paste0(
    m$desc, " of the CATE intervals, confidence-interval study, ", o$label,
    ", scenarios 1--4, by the subsampling ratio of the half-sample bootstrap: ",
    cells,
    "Each model's interval is its 95\\% half-sample bootstrap band, ",
    "simultaneous over units. ", m$defn,
    "Causal forest (inbuilt CI): the pointwise 95\\% normal interval from the ",
    "causal forest's own variance estimate, which does not depend on the ",
    "subsampling ratio, so it is printed once."
  )
}

sf_caption <- function(o, runs) {
  paste0(
    "Data-driven subsampling ratio, DR-RandomForest, ", o$label,
    ", correlated covariates, scenarios 1--4: mean (Monte Carlo SE) over ",
    runs, " runs per cell. Each run picks the ratio whose 90\\% band covers ",
    "the run's own estimate closest to 90\\% of the time, and the band is ",
    "then built at that ratio. Plug-in coverage is that calibration's own ",
    "coverage at the pick; marginal and simultaneous coverage are of the ",
    "true CATE (nominal 0.90). $\\Delta$: the per-run difference, ",
    "$\\rho = 0.5$ minus $\\rho = 0$; a run shares its random draws at the two ",
    "values of $\\rho$, so its SE is that of the paired difference."
  )
}

# ---- the CI_sf sweep ----------------------------------------------------------

sweep_design <- merge(DESIGN, data.frame(CI_sf = CI_SF))

for (o in sets) {
  metrics <- read_complete(
    file.path(res_root, o$dir, paste0(o$prefix, "_corr_ci_metrics.RDS")),
    sweep_design, c("rho", "scenario", "n", "CI_sf"), o$out
  )
  if (is.null(metrics)) next

  metrics <- filter(metrics, !grepl("_grid$", model))
  keys <- c("scenario", "n", "CI_sf", "model")
  metrics <- if (o$diff) {
    paired_rho_diff(metrics, ci_metrics$col, c(keys, "run"))
  } else {
    filter(metrics, rho == o$rho)
  }

  metrics_summary <- summarise_metrics(
    metrics, keys,
    cols = setNames(ci_metrics$col, ci_metrics$stem),
    binomial = if (o$diff) character() else ci_metrics$stem[ci_metrics$binomial],
    count_na = character()
  ) %>%
    apply_labels(SS_SCENARIO_LABELS) %>%
    mutate(n = factor(n, levels = c(500, 1000)))

  runs <- runs_per_cell(metrics, keys)

  for (i in seq_len(nrow(ci_metrics))) {
    m <- ci_metrics[i, ]
    digits <- if (is.na(m$digits)) o$digits_len else m$digits
    name <- paste0(o$out, "_ci_", m$file)
    tex <- ci_sweep_table(
      metrics_summary, m$stem, digits,
      span_models = SPAN_MODELS,
      caption.short = paste0(m$desc, ", CATE intervals, ", o$label),
      caption = sweep_caption(o, m, runs),
      label = name
    )
    out_file <- file.path(tab_path, paste0(name, ".tex"))
    writeLines(tex, out_file)
    message("wrote ", out_file)
  }
}

# ---- the data-driven CI_sf ----------------------------------------------------

# column stem, source column, header, digits (NA: the outcome's length digits)
sf_cols <- tibble::tribble(
  ~stem,        ~col,                    ~header,                                          ~digits, ~binomial,
  "pick",       "optimal_sf",            "\\shortstack[r]{Ratio\\\\picked}",               2,       FALSE,
  "plugin_cov", "plugin_coverage",       "\\shortstack[r]{Plug-in\\\\coverage}",           3,       FALSE,
  "marg_cov",   "marginal_coverage",     "\\shortstack[r]{Marginal\\\\coverage}",          3,       FALSE,
  "simul_cov",  "simultaneous_coverage", "\\shortstack[r]{Simultaneous\\\\coverage}",      2,       TRUE,
  "ci_len",     "mean_ci_length",        "\\shortstack[r]{Mean\\\\length}",                NA,      FALSE
)

for (o in outcomes[1]) {  # the binary optimal_sf study was dropped
  name <- paste0(o$prefix, "_corr_ci_sf")
  metrics <- read_complete(
    file.path(res_root, o$dir, "sf_calibration", paste0(o$prefix, "_corr_ci_sf_metrics.RDS")),
    DESIGN, c("rho", "scenario", "n"), name
  )
  if (is.null(metrics)) next

  # its single arm, under the name its metrics script files it as
  metrics <- filter(metrics, model == "dr_random_forest")
  keys <- c("scenario", "n")
  summ <- function(df, binomial) {
    summarise_metrics(df, keys, cols = setNames(sf_cols$col, sf_cols$stem),
                      binomial = binomial, count_na = character())
  }
  binom <- sf_cols$stem[sf_cols$binomial]
  sf_summary <- bind_rows(
    summ(filter(metrics, rho == 0), binom) %>% mutate(rho = "0"),
    summ(filter(metrics, rho == 0.5), binom) %>% mutate(rho = "0.5"),
    summ(paired_rho_diff(metrics, sf_cols$col, c(keys, "run")), character()) %>%
      mutate(rho = "$\\Delta$")
  ) %>%
    apply_labels(SS_SCENARIO_LABELS) %>%
    mutate(n = factor(n, levels = c(500, 1000)),
           rho = factor(rho, levels = c("0", "0.5", "$\\Delta$"))) %>%
    arrange(scenario, n, rho)

  cells <- lapply(seq_len(nrow(sf_cols)), function(i) {
    digits <- if (is.na(sf_cols$digits[i])) o$digits_len else sf_cols$digits[i]
    fmt_mcse(sf_summary[[paste0("mean_", sf_cols$stem[i])]],
             sf_summary[[paste0("mcse_", sf_cols$stem[i])]], digits)
  })
  names(cells) <- sf_cols$stem

  tex <- grouped_longtable(
    bind_cols(tibble(n = as.character(sf_summary$n),
                     rho = as.character(sf_summary$rho)),
              as_tibble(cells)),
    col_names = c("$n$", "$\\rho$", sf_cols$header),
    group = sf_summary$scenario,
    block = sf_summary$n,
    caption.short = paste0("Data-driven subsampling ratio, ", o$label),
    caption = sf_caption(o, runs_per_cell(metrics, c("rho", keys))),
    label = name,
    # five columns fit portrait at 10pt, as ss_tables.R's four
    landscape = FALSE,
    font_size = 10,
    tabcolsep = NULL
  )
  out_file <- file.path(tab_path, paste0(name, ".tex"))
  writeLines(tex, out_file)
  message("wrote ", out_file)

  # the main-text float: the rho = 0.5 rows, with the pick as median (IQR).
  # type = 1 quantiles stay on the ratio grid, where median() would average
  # the two middle picks of 100 runs
  main <- filter(metrics, rho == 0.5)
  picks <- main %>%
    group_by(scenario, n) %>%
    summarise(q1 = quantile(optimal_sf, 0.25, type = 1),
              med = quantile(optimal_sf, 0.5, type = 1),
              q3 = quantile(optimal_sf, 0.75, type = 1),
              .groups = "drop") %>%
    mutate(pick = sprintf("%.2f (%.2f--%.2f)", med, q1, q3)) %>%
    apply_labels(SS_SCENARIO_LABELS) %>%
    mutate(n = factor(n, levels = c(500, 1000)))

  main_rows <- sf_summary$rho == "0.5"
  main_stems <- setdiff(sf_cols$stem, "pick")
  body <- sf_summary[main_rows, c("scenario", "n")] %>%
    left_join(select(picks, scenario, n, pick), by = c("scenario", "n"))
  first_sc <- !duplicated(body$scenario)
  tab <- cbind(
    ifelse(seq_len(nrow(body)) == 1, tools::toTitleCase(sub(" outcome$", "", o$label)), ""),
    ifelse(first_sc, as.character(body$scenario), ""),
    as.character(body$n),
    body$pick,
    sapply(cells[main_stems], `[`, main_rows)
  )
  header <- c("Outcome", "Scenario", "$n$", "\\shortstack[r]{Ratio, median\\\\(IQR)}",
              sf_cols$header[match(main_stems, sf_cols$stem)])
  rows <- paste0(ifelse(first_sc & seq_along(first_sc) > 1, "\\addlinespace\n", ""),
                 apply(tab, 1, paste, collapse = " & "), " \\\\")

  main_caption <- paste0(
    "Data-driven choice of the subsampling ratio for the half-sample bootstrap ",
    "bands of the \\gls{dr}-random forest, targeting 90\\% coverage, correlated ",
    "covariates ($\\rho = 0.5$) \\par \\small Ratio: the median (interquartile ",
    "range) over runs of the ratio each run's calibration picked. Plug-in ",
    "coverage: the calibration's mean coverage of $\\hat{\\tau}(X)$ at the ",
    "selected ratio. Coverage of the true \\gls{cate} and length are for the ",
    "final band ($B = 200$). Other cells: mean (Monte Carlo SE) over ",
    runs_per_cell(main, keys), " runs per row."
  )
  # tabcolsep at grouped_longtable()'s 3pt, local to the float: at LaTeX's 6pt
  # the eight columns overrun a 16cm text width by about 1cm
  tex <- paste0(
    "\\begin{table}[ht]\n\\centering\n\\small\n\\setlength{\\tabcolsep}{3pt}\n",
    "\\caption[Data-driven choice of the subsampling ratio]{", main_caption, "}\n",
    "\\label{tab:optimal_sf}\n",
    "\\begin{tabular}{@{}lll", strrep("r", length(header) - 3), "@{}}\n",
    "\\toprule\n",
    paste(header, collapse = " & "), " \\\\ \\midrule\n",
    paste(rows, collapse = "\n"), "\n",
    "\\bottomrule\n\\end{tabular}\n\\end{table}\n"
  )
  out_file <- file.path(tab_path, paste0(name, "_main.tex"))
  writeLines(tex, out_file)
  message("wrote ", out_file)
}
