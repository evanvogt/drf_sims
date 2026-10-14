##########
# title: LaTeX tables for the thesis appendix - model evaluation
##########
# The model-evaluation study's strategies tree (model_evaluation/), the
# appendix counterpart of thesis_figures/model_eval.R: the same 24 proxies (DR
# and IF, estimated or fixed propensity, three fold regimes, two machines),
# the same per-run selection and the same references, from
# model_evaluation/me_selection.R. Labels follow the figures' legends.
#
# Measures, per (scenario, n, run, proxy), each "mean (MCSE)" over runs:
#   - top-1 selection accuracy, whether the pick is the truly best candidate
#     (binomial MCSE, as in the figure);
#   - the pick's true rank among the 9 candidates, 1 = best;
#   - regret, the pick's true PEHE minus the lowest in the run.
#
# Tables:
#   - one per measure, rows scenario / n / score type, columns fold regime x
#     machine, with a random-pick and a best-fixed-choice row under each n
#     (me_proxy_table(), R/tables.R). The best fixed choice is the candidate
#     with the lowest mean regret, used for all three measures;
#   - one design-axis table, each axis level averaged over the other axes
#     within run, then over runs, as the axis figure does (me_axis_runs()).
#     Top-1 averaged within a run is a proportion, so there it takes the
#     general MCSE, not the binomial one;
#   - one fixed-choice table, rows scenario / n / candidate: each candidate's
#     true PEHE and the three measures if it were always deployed
#     (me_candidates()), the best fixed choice in bold. Its bold rows are the
#     proxy tables' best-fixed-choice rows;
#   - one completed-runs table, the proxy tables' layout with each cell the
#     number of runs that proxy was scored in, and an "All runs" row (the runs
#     with a true PEHE) under each n.
#
# Run counts are read off the data, not the design. A proxy whose scoring
# nuisance failed in a run is NA for every candidate there, and me_selection()
# drops that (run, proxy), so a failure shows only as a short count. Before
# the scoring propensity was trimmed (DEVLOG bug S) that was AutoML's
# estimated-pi proxies, mostly holdout; on metrics computed since, every count
# should be full, and this is the check that it is. Each caption gives its own
# table's counts, and the proxy tables name any short proxies and point to the
# completed-runs table.
#
# All six are portrait, within the thesis's 17cm text width. Regret's
# decimals are read off its MCSEs (mcse_digits() below): its scale changes
# with the scenario and n, and the results were not to hand when this was
# written.
#
# Writes to ../results/thesis_tables/:
#   me_top1.tex, me_rank.tex, me_regret.tex, me_axis.tex, me_fixed.tex,
#   me_runs.tex
# \input{} them into a document with
#   \usepackage{booktabs, longtable, pdflscape, array}
#
# Optional arguments, for a smoke run on mock metrics:
#   Rscript me_tables.R [metrics.RDS] [table dir]

library(here)
source(here("R", "figures.R"))
source(here("R", "tables.R"))
source(here("model_evaluation", "me_selection.R"))

# paths
args <- commandArgs(trailingOnly = TRUE)
path <- here()
metrics_file <- if (length(args) >= 1) {
  args[1]
} else {
  file.path(dirname(path), "results", "correlated",
            "model_evaluation_strategies", "me_strat_metrics.RDS")
}
tab_path <- if (length(args) >= 2) {
  args[2]
} else {
  file.path(dirname(path), "results", "thesis_tables")
}
dir.create(tab_path, showWarnings = FALSE, recursive = TRUE)

# display labels, as thesis_figures/model_eval.R's but in LaTeX (the tables
# are written with escape = FALSE). The order of each vector is the row or
# column order.
SCORE_LABELS <- c(dr = "DR", dr05 = "DR-$\\pi$-fixed", infl = "IF",
                  infl05 = "IF-$\\pi$-fixed")
FOLD_LABELS <- c(cv_shared = "CV", holdout = "Holdout", whole = "Whole")
PIPELINE_LABELS <- c(xgb = "XGBoost", automl = "AutoML")

AXES <- list(
  Score = list(col = "score_family", labels = c(dr = "DR", infl = "IF")),
  Propensity = list(col = "pi_regime",
                    labels = c(estimated = "Estimated", fixed = "Fixed at 0.5")),
  `Fold regime` = list(col = "fold_type", labels = FOLD_LABELS),
  Machine = list(col = "pipeline", labels = PIPELINE_LABELS)
)

MEASURE_COLS <- c(hit = "hit", rank = "pick_rank", regret = "regret")

#' Decimals for two significant figures of the typical (median) MCSE, from 2
#' to 4: a fifth would not fit the proxy tables' six columns in portrait
mcse_digits <- function(mcse) {
  m <- median(mcse[mcse > 0], na.rm = TRUE)
  if (is.na(m)) return(2)
  min(4, max(2, 1 - floor(log10(m))))
}

# ---- data --------------------------------------------------------------------

metrics <- me_read_metrics(metrics_file, SS_SCENARIO_LABELS)
selection <- me_selection(metrics)
candidates <- me_candidates(metrics)

# ---- run counts --------------------------------------------------------------

#' "358", or "289--350" when the counts differ
count_range <- function(x) {
  r <- range(x)
  if (r[1] == r[2]) as.character(r[1]) else paste0(r[1], "--", r[2])
}

# a failed candidate would shrink the "out of 9" the captions state, so say so
short_runs <- candidates %>%
  count(scenario, n, run, name = "candidates") %>%
  filter(candidates < ME_N_CANDIDATES)
if (nrow(short_runs) > 0) {
  warning(nrow(short_runs), " run(s) have fewer than ", ME_N_CANDIDATES,
          " candidates with a true PEHE; the captions' \"of ",
          ME_N_CANDIDATES, "\" is then an upper bound", call. = FALSE)
}

# the runs with a true PEHE, per (scenario, n): every run the reference rows,
# the axis table and the fixed-choice table can use
total_runs <- candidates %>%
  distinct(scenario, n, run) %>%
  count(scenario, n, name = "total")

# over the full proxy grid, so a proxy that failed in every run counts 0
# rather than vanishing
proxy_runs <- tidyr::expand_grid(
  total_runs,
  score_type = ME_SCORE_TYPES, fold_type = ME_FOLD_TYPES,
  pipeline = ME_PIPELINES
) %>%
  left_join(count(selection, scenario, n, score_type, fold_type, pipeline,
                  name = "runs"),
            by = c("scenario", "n", "score_type", "fold_type", "pipeline")) %>%
  mutate(runs = coalesce(runs, 0L))

short_cells <- proxy_runs %>% filter(runs < total)
short_proxies <- short_cells %>%
  distinct(fold_type, pipeline) %>%
  mutate(label = paste0(FOLD_LABELS[fold_type], "/",
                        PIPELINE_LABELS[pipeline]))
short_proxies <- paste(short_proxies$label, collapse = " and ")

if (nrow(short_cells) > 0) {
  message("proxy cells with failed runs (", nrow(short_cells), " of ",
          nrow(proxy_runs), "):")
  print(as.data.frame(short_cells), row.names = FALSE)
} else {
  message("no proxy cell has failed runs")
}

proxy_runs_text <- paste0(count_range(total_runs$total), " runs per cell")
if (nrow(short_cells) > 0) {
  proxy_runs_text <- paste0(
    proxy_runs_text, ", fewer for the ", short_proxies, " proxies (",
    count_range(short_cells$runs), "; Table~\\ref{tab:me_runs})"
  )
}

# ---- summaries ---------------------------------------------------------------

proxy_summary <- selection %>%
  mutate(
    metric = label_factor(score_type, SCORE_LABELS),
    fold = label_factor(fold_type, FOLD_LABELS),
    machine = label_factor(pipeline, PIPELINE_LABELS)
  ) %>%
  summarise_metrics(
    c("scenario", "n", "metric", "fold", "machine"),
    cols = MEASURE_COLS,
    binomial = "hit",
    count_na = character()
  )

# the random pick's hit is a run's mean over its candidates, not 0/1, so only
# the best fixed choice's takes the binomial MCSE
ref_runs <- me_references(metrics)
ref_summary <- bind_rows(
  ref_runs %>%
    filter(ref == "random") %>%
    summarise_metrics(c("scenario", "n", "ref", "model"),
                      cols = MEASURE_COLS, count_na = character()),
  ref_runs %>%
    filter(ref == "best_fixed") %>%
    summarise_metrics(c("scenario", "n", "ref", "model"),
                      cols = MEASURE_COLS, binomial = "hit",
                      count_na = character())
) %>%
  mutate(label = ifelse(ref == "random", "Random pick",
                        paste0("Best fixed (", model, ")")))

axis_runs <- me_axis_runs(selection, AXES, unname(MEASURE_COLS))
axis_summary <- axis_runs %>%
  summarise_metrics(c("scenario", "n", "axis", "level"),
                    cols = MEASURE_COLS, count_na = character())

regret_digits <- mcse_digits(c(proxy_summary$mcse_regret,
                               axis_summary$mcse_regret))
message("regret to ", regret_digits, " decimals")

# ---- captions ----------------------------------------------------------------

proxy_key <- paste0(
  "Rows: the proxy's score, DR risk (DR) or the influence-corrected PEHE ",
  "(IF), with the scoring nuisance's estimated propensity or with it fixed ",
  "at the true 0.5 ($\\pi$-fixed). Columns: the rows the scoring nuisance is ",
  "fit on, the candidate's training folds, predicting its held-out fold (CV), ",
  "the held-out fold alone (Holdout; at $n = 250$ adjacent pairs of folds are ",
  "pooled), or all $n$ rows (Whole); and its machine. Random pick: a ",
  "candidate drawn uniformly at random. Best fixed: always deploying the ",
  "named candidate, the one with the lowest mean regret in that scenario and ",
  "$n$, chosen with hindsight."
)

MEASURES <- list(
  hit = list(
    file = "me_top1", digits = 2, short = "Top-1 selection accuracy",
    text = paste0(
      "the proportion of runs in which the proxy's pick is the candidate ",
      "with the lowest true PEHE, out of 9 (a random pick's is 1/9)"
    )
  ),
  rank = list(
    file = "me_rank", digits = 2, short = "True rank of the selected candidate",
    text = paste0(
      "the true rank of the proxy's pick among the 9 candidates, 1 being the ",
      "lowest true PEHE, tied candidates sharing the better rank (a random ",
      "pick's is 5)"
    )
  ),
  regret = list(
    file = "me_regret", digits = regret_digits, short = "Regret",
    text = paste0(
      "the true PEHE of the proxy's pick minus the lowest true PEHE among the ",
      "candidates in the run"
    )
  )
)

caption_text <- function(m) {
  paste0(
    "Model evaluation, ", tolower(substr(m$short, 1, 1)), substring(m$short, 2),
    " of the selection proxies: ", m$text, ". Mean (Monte Carlo SE) over ",
    proxy_runs_text, ". ", proxy_key
  )
}

# ---- tables ------------------------------------------------------------------

# 8pt with 3pt column padding: the six "mean (MCSE)" columns then fit 17cm in
# portrait at up to 4 decimals
for (stem in names(MEASURES)) {
  m <- MEASURES[[stem]]
  tex <- me_proxy_table(
    proxy_summary, ref_summary, stem, m$digits,
    caption.short = paste0("Model evaluation: ", tolower(m$short)),
    caption = caption_text(m),
    label = m$file,
    landscape = FALSE,
    font_size = 8,
    tabcolsep = "3pt"
  )
  out_file <- file.path(tab_path, paste0(m$file, ".tex"))
  writeLines(tex, out_file)
  message("wrote ", out_file)
}

# two-line headers keep the columns as narrow as their cells
axis_cols <- tibble::tribble(
  ~stem,    ~header,                                   ~digits,
  "hit",    "\\shortstack[r]{Top-1\\\\accuracy}",      2,
  "rank",   "\\shortstack[r]{True rank\\\\of pick}",   2,
  "regret", "Regret",                                  regret_digits
)

# a run keeps an axis level while any proxy sharing it was scored, so a failed
# proxy changes what the level averages over in that run, not the run count
axis_missing <- if (nrow(short_cells) > 0) {
  paste0(
    " The ", short_proxies, " proxies are missing from some runs ",
    "(Table~\\ref{tab:me_runs}); in those runs, the levels they share are ",
    "averaged over the proxies that remain."
  )
} else {
  ""
}

tex <- me_axis_table(
  axis_summary, axis_cols,
  caption.short = "Model evaluation: selection by design axis",
  caption = paste0(
    "Model evaluation, selection by design axis: top-1 accuracy, true rank ",
    "of the pick (1 = best of 9) and regret, each level of an axis averaged ",
    "over the proxies that share it within each run, then over runs. Mean ",
    "(Monte Carlo SE) over ",
    runs_per_cell(axis_runs, c("scenario", "n", "axis", "level")),
    " runs per cell; the SE counts runs. ",
    "Averaged within a run, top-1 accuracy is a proportion rather than 0/1, ",
    "so its SE is the general one. Fold regime: the rows the scoring ",
    "nuisance is fit on, as in Table~\\ref{tab:me_top1}.", axis_missing
  ),
  label = "me_axis",
  landscape = FALSE,
  font_size = 10,
  tabcolsep = NULL
)
out_file <- file.path(tab_path, "me_axis.tex")
writeLines(tex, out_file)
message("wrote ", out_file)

# ---- the candidates as fixed choices -----------------------------------------

# the families in model_evaluation/ADEMP.md's order: forests, glmnet, SL
CANDIDATE_LEVELS <- c("rf1", "rf2", "rf3", "net1", "net2", "net3",
                      "SL1", "SL2", "SL3")

best_fixed <- ref_summary %>%
  filter(ref == "best_fixed") %>%
  select(scenario, n, best = model)

# the best fixed choice's name is bolded in place, so each candidate has a
# plain and a bold level, kept next to each other in the row order
candidate_summary <- candidates %>%
  summarise_metrics(
    c("scenario", "n", "model"),
    cols = c(pehe = "true_pehe", MEASURE_COLS),
    binomial = "hit",
    count_na = character()
  ) %>%
  left_join(best_fixed, by = c("scenario", "n")) %>%
  mutate(model = as.character(model))
ids <- union(intersect(CANDIDATE_LEVELS, candidate_summary$model),
             sort(unique(candidate_summary$model)))
candidate_summary <- candidate_summary %>%
  mutate(model = factor(
    ifelse(model == best, paste0("\\textbf{", model, "}"), model),
    levels = as.vector(rbind(ids, paste0("\\textbf{", ids, "}")))
  ))

# the axis table's columns, after the true PEHE; a candidate is not a pick,
# so its rank is just its true rank
fixed_cols <- bind_rows(
  tibble::tibble(stem = "pehe", header = "\\shortstack[r]{True\\\\PEHE}",
                 digits = mcse_digits(candidate_summary$mcse_pehe)),
  axis_cols
) %>%
  mutate(header = ifelse(stem == "rank", "\\shortstack[r]{True\\\\rank}",
                         header))

tex <- ss_latex_table(
  candidate_summary, fixed_cols,
  caption.short = "Model evaluation: the candidates as fixed choices",
  caption = paste0(
    "Model evaluation, the candidates as fixed choices: each candidate's ",
    "true PEHE, and its top-1 accuracy (the proportion of runs in which it is ",
    "the candidate with the lowest true PEHE), true rank (1 = best of 9) and ",
    "regret if it were always deployed. Mean (Monte Carlo SE) over ",
    runs_per_cell(candidates, c("scenario", "n", "model")),
    " runs per cell. In bold, the best fixed choice: the candidate with the ",
    "lowest mean regret in that scenario and $n$, chosen with hindsight, and ",
    "the Best fixed row of Tables~\\ref{tab:me_top1}--\\ref{tab:me_regret}."
  ),
  label = "me_fixed",
  landscape = FALSE,
  font_size = 10,
  tabcolsep = NULL
)
out_file <- file.path(tab_path, "me_fixed.tex")
writeLines(tex, out_file)
message("wrote ", out_file)

# ---- completed runs ----------------------------------------------------------

# the proxy tables' layout, so each count sits where the cell it is behind
# does. A count has no MCSE, which fmt_mcse() prints as the number alone.
runs_summary <- proxy_runs %>%
  mutate(
    metric = label_factor(score_type, SCORE_LABELS),
    fold = label_factor(fold_type, FOLD_LABELS),
    machine = label_factor(pipeline, PIPELINE_LABELS),
    mean_runs = runs,
    mcse_runs = NA_real_
  )
runs_refs <- total_runs %>%
  mutate(ref = factor("total"), label = "All runs", mean_runs = total,
         mcse_runs = NA_real_)

tex <- me_proxy_table(
  runs_summary, runs_refs, "runs", 0,
  caption.short = "Model evaluation: completed runs per proxy",
  caption = paste0(
    "Model evaluation, completed runs: the number of runs in which each ",
    "selection proxy could be scored, the runs behind each cell of ",
    "Tables~\\ref{tab:me_top1}--\\ref{tab:me_regret}. A proxy is missing from ",
    "a run when its scoring nuisance returned missing values there, which ",
    "leaves the proxy undefined for every candidate. All runs: the runs with ",
    "a true PEHE for the candidates, the runs behind the Random pick and Best ",
    "fixed rows and Tables~\\ref{tab:me_axis} and \\ref{tab:me_fixed}. Rows ",
    "and columns as in Table~\\ref{tab:me_top1}."
  ),
  label = "me_runs",
  landscape = FALSE,
  font_size = 10,
  tabcolsep = NULL
)
out_file <- file.path(tab_path, "me_runs.tex")
writeLines(tex, out_file)
message("wrote ", out_file)
