##########
# title: LaTeX tables for the thesis appendix - model evaluation
##########
# Appendix tables for the model-evaluation strategies tree, the counterpart of
# thesis_figures/model_eval.R: the same 24 proxies, per-run selection and
# references (model_evaluation/me_selection.R). Cells are "mean (MCSE)" over
# runs of top-1 accuracy, the pick's true rank, or regret.
#
# - The best fixed choice is the candidate with the lowest mean regret, used for
#   all three measures.
# - Design-axis levels are averaged within run, then over runs, so top-1 there
#   is a proportion and takes the general MCSE, not the binomial one.
# - Run counts are read off the data: a proxy whose nuisance failed in a run
#   drops out of that run, so a short count flags a failure (me_runs.tex).
#
# Writes to ../results/thesis_tables/: me_top1, me_rank, me_regret, me_axis,
# me_fixed, me_runs, me_if_population, me_if_selection (.tex). me_if_population
# needs model_evaluation/me_infl_investigation.R's population.RDS and is skipped
# without it. \input{} with \usepackage{booktabs, longtable, pdflscape, array}.
#
#   Rscript me_tables.R [metrics.RDS] [table dir] [population.RDS]   # mock smoke run

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
population_file <- if (length(args) >= 3) {
  args[3]
} else {
  file.path(dirname(path), "results", "correlated", "model_evaluation",
            "investigation", "population.RDS")
}

# display labels, as thesis_figures/model_eval.R's but in LaTeX (the tables
# are written with escape = FALSE). The order of each vector is the row or
# column order.
SCORE_LABELS <- c(dr = "DR", dr05 = "DR-$\\pi$-fixed", inflt = "IF",
                  inflt05 = "IF-$\\pi$-fixed")
FOLD_LABELS <- c(cv_shared = "CV", holdout = "Holdout", whole = "Whole")
PIPELINE_LABELS <- c(xgb = "XGBoost", automl = "AutoML")

AXES <- list(
  Score = list(col = "score_family", labels = c(dr = "DR", inflt = "IF")),
  `Propensity score` = list(col = "pi_regime",
                            labels = c(estimated = "estimated",
                                       fixed = "fixed at 0.5")),
  `Data splitting` = list(col = "fold_type", labels = FOLD_LABELS),
  `Surrogate learner` = list(col = "pipeline", labels = PIPELINE_LABELS)
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

# the IF-selection table's extra columns: the phi-surrogate IF score in its
# arm, and the MSE against the true mu1. Every other table ignores them.
IF_SEL_ARM <- "cv_shared_xgb"
IF_SEL_SCORES <- c(dr = "dr05", infl = "infl05", inflt = "inflt05")
IF_SEL_COLS <- setNames(paste0(IF_SEL_SCORES, "_", IF_SEL_ARM),
                        names(IF_SEL_SCORES))

metrics <- me_read_metrics(metrics_file, SS_SCENARIO_LABELS,
                           extra_cols = c(IF_SEL_COLS[["infl"]], "mu1_mse"))
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

proxy_runs_text <- paste0(count_range(total_runs$total),
                          " simulated trials per cell")
if (nrow(short_cells) > 0) {
  proxy_runs_text <- paste0(
    proxy_runs_text, ", fewer for the ", short_proxies, " metrics (",
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
  mutate(label = ifelse(ref == "random", "Random selection",
                        paste0("Best fixed choice (", model, ")")))

axis_runs <- me_axis_runs(selection, AXES, unname(MEASURE_COLS))
axis_summary <- axis_runs %>%
  summarise_metrics(c("scenario", "n", "axis", "level"),
                    cols = MEASURE_COLS, count_na = character())

regret_digits <- mcse_digits(c(proxy_summary$mcse_regret,
                               axis_summary$mcse_regret))
message("regret to ", regret_digits, " decimals")

# ---- captions ----------------------------------------------------------------

proxy_key <- paste0(
  "Rows: score (DR: the DR score; IF: the influence-function-corrected ",
  "\\gls{pehe} with the T-learner surrogate $\\tilde{\\mu}_1 - ",
  "\\tilde{\\mu}_0$), with the propensity score estimated by the surrogate ",
  "pipeline and truncated to $[0.05, 0.95]$, or fixed at the randomisation ",
  "probability of 0.5 ($\\pi$-fixed). Columns: data splitting for the ",
  "surrogate's nuisance functions (CV: fitted on the candidate's training ",
  "folds and predicted in its held-out fold; Holdout: fitted and predicted ",
  "within the held-out fold, with adjacent pairs of folds pooled at ",
  "$n = 250$; Whole: fitted and predicted in all $n$ participants), and ",
  "surrogate nuisance learner (XGBoost or AutoML). Random selection: a ",
  "candidate chosen uniformly at random. Best fixed choice: always selecting ",
  "the named candidate, the one with the lowest mean regret in that scenario ",
  "and $n$, identified with hindsight. Mean (Monte Carlo SE) over ",
  proxy_runs_text, "."
)

MEASURES <- list(
  hit = list(
    file = "me_top1", digits = 2, short = "Top-1 selection accuracy",
    text = paste0(
      "Top-1 selection accuracy of each metric: the proportion of simulated ",
      "trials in which the selected candidate had the lowest true \\gls{pehe} ",
      "of the nine (1/9 for random selection)."
    )
  ),
  rank = list(
    file = "me_rank", digits = 2, short = "True rank of the selected candidate",
    text = paste0(
      "True rank of the candidate selected by each metric, from 1 (lowest ",
      "true \\gls{pehe}) to 9, with tied candidates sharing the better rank ",
      "(5 for random selection)."
    )
  ),
  regret = list(
    file = "me_regret", digits = regret_digits, short = "Regret",
    text = paste0(
      "Regret of each metric: the true \\gls{pehe} of the selected candidate ",
      "minus the lowest true \\gls{pehe} among the nine candidates in the ",
      "same simulated trial."
    )
  )
)

caption_text <- function(m) paste(m$text, proxy_key)

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
  "rank",   "\\shortstack[r]{True rank\\\\of selected}", 2,
  "regret", "Regret",                                  regret_digits
)

# a run keeps an axis level while any proxy sharing it was scored, so a failed
# proxy changes what the level averages over in that run, not the run count
axis_missing <- if (nrow(short_cells) > 0) {
  paste0(
    " The ", short_proxies, " metrics are missing from some simulated ",
    "trials (Table~\\ref{tab:me_runs}); in those trials, the levels they ",
    "share are averaged over the metrics that remain."
  )
} else {
  ""
}

tex <- me_axis_table(
  axis_summary, axis_cols,
  caption.short = "Model evaluation: selection by design factor",
  caption = paste0(
    "Selection performance by design factor: top-1 accuracy, true rank of ",
    "the selected candidate (1 = best of 9) and regret, for each level of ",
    "one factor, averaged first over the metrics that share that level ",
    "within each simulated trial and then over simulated trials. Because it ",
    "is averaged within a trial, top-1 accuracy is a proportion rather than ",
    "0 or 1. Factors as in Table~\\ref{tab:me_design}. Mean (Monte Carlo SE) ",
    "over ", runs_per_cell(axis_runs, c("scenario", "n", "axis", "level")),
    " simulated trials per cell.", axis_missing
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
    "Each candidate as a fixed choice: its true \\gls{pehe}, its top-1 ",
    "accuracy (the proportion of simulated trials in which it had the lowest ",
    "true \\gls{pehe}), its true rank (1 = best of 9) and its regret if it ",
    "were always selected. Mean (Monte Carlo SE) over ",
    runs_per_cell(candidates, c("scenario", "n", "model")),
    " simulated trials per cell. In bold: the best fixed choice, the ",
    "candidate with the lowest mean regret in that scenario and $n$, ",
    "identified with hindsight, which is the Best fixed choice row of ",
    "Tables~\\ref{tab:me_top1}--\\ref{tab:me_regret}."
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
  mutate(ref = factor("total"), label = "All simulated trials",
         mean_runs = total, mcse_runs = NA_real_)

tex <- me_proxy_table(
  runs_summary, runs_refs, "runs", 0,
  caption.short = "Model evaluation: completed simulated trials per metric",
  caption = paste0(
    "Completed simulated trials: the number of simulated trials in which ",
    "each metric could be computed, behind each cell of ",
    "Tables~\\ref{tab:me_top1}--\\ref{tab:me_regret}. A metric is missing ",
    "from a simulated trial when the surrogate's nuisance functions returned ",
    "missing values there, which leaves the metric undefined for every ",
    "candidate. All simulated trials: those with a true \\gls{pehe} for the ",
    "candidates, behind the Random selection and Best fixed choice rows and ",
    "Tables~\\ref{tab:me_axis} and \\ref{tab:me_fixed}. Rows and columns as ",
    "in Table~\\ref{tab:me_top1}."
  ),
  label = "me_runs",
  landscape = FALSE,
  font_size = 10,
  tabcolsep = NULL
)
out_file <- file.path(tab_path, "me_runs.tex")
writeLines(tex, out_file)
message("wrote ", out_file)

# ---- the IF score: population check ------------------------------------------

# the constructed candidates of me_infl_investigation.R's investigate_scenario(),
# in the order that breaks a tie in true PEHE
POP_CANDIDATE_LABELS <- c(
  truth = "$\\tau$",
  half = "$\\bar\\tau + (\\tau - \\bar\\tau)/2$",
  ate = "$\\bar\\tau$",
  mu1 = "$\\mu_1$",
  tau_plus_X2 = "$\\tau + X_2$"
)

# every column a loss. With the oracle nuisances the T-learner surrogate
# mu1 - mu0 is the true tau, so infl_tau is the IF score with the T-learner
# surrogate.
pop_cols <- tibble::tribble(
  ~stem,       ~header,
  "true_pehe", "\\shortstack[r]{True\\\\PEHE}",
  "dr",        "DR",
  "infl_tau",  "\\shortstack[r]{IF,\\\\true $\\tau$}",
  "infl_phi",  "\\shortstack[r]{IF,\\\\pseudo-outcome}",
  "mu1_mse",   "\\shortstack[r]{MSE\\\\vs $\\mu_1$}"
)
POP_DIGITS <- 3

if (file.exists(population_file)) {
  population <- readRDS(population_file)
  absent <- setdiff(c("n_pop", pop_cols$stem), names(population))
  if (length(absent) > 0) {
    stop(population_file, " has no ", paste(absent, collapse = ", "),
         ". Rerun `Rscript me_infl_investigation.R` in model_evaluation/.",
         call. = FALSE)
  }

  # ranked on the printed values, so candidates that differ only past the
  # printed decimals (tau, its mean and the half-way candidate are the same
  # vector in the null scenario, up to rounding) share a rank, the better one
  pop <- population %>%
    filter(as.character(scenario) %in% names(SS_SCENARIO_LABELS)) %>%
    mutate(
      scenario = label_factor(scenario, SS_SCENARIO_LABELS),
      candidate = label_factor(candidate, POP_CANDIDATE_LABELS)
    ) %>%
    group_by(scenario) %>%
    mutate(across(
      all_of(pop_cols$stem),
      ~ rank(round(.x, POP_DIGITS), ties.method = "min"),
      .names = "rank_{.col}"
    )) %>%
    ungroup() %>%
    arrange(scenario, rank_true_pehe, candidate)

  pop_cells <- lapply(pop_cols$stem, function(s) {
    paste0(fmt_mcse(pop[[s]], NA_real_, POP_DIGITS),
           " (", pop[[paste0("rank_", s)]], ")")
  })
  names(pop_cells) <- pop_cols$stem

  n_pop <- unique(pop$n_pop)
  if (length(n_pop) > 1) {
    stop("population.RDS mixes sample sizes: ", paste(n_pop, collapse = ", "),
         call. = FALSE)
  }
  n_pop_text <- paste0("$N = ",
                       formatC(n_pop, format = "d", big.mark = "{,}"), "$")

  tex <- grouped_longtable(
    bind_cols(tibble::tibble(candidate = as.character(pop$candidate)),
              tibble::as_tibble(pop_cells)),
    col_names = c("Candidate", pop_cols$header),
    group = pop$scenario,
    block = NULL,
    align = c("l", rep("r", nrow(pop_cols))),
    caption.short = paste0("Model evaluation: the IF score with the true ",
                           "nuisance functions"),
    caption = paste0(
      "The IF score with the true nuisance functions: five constructed ",
      "candidates scored in one sample of ", n_pop_text, " from each ",
      "scenario, using the true $\\mu_0$, $\\mu_1$ and $\\pi = 0.5$, so that ",
      "any departure from the true \\gls{pehe} is due to the score itself ",
      "rather than to error in the nuisance functions. Each column is a ",
      "score; in brackets, its rank among the five candidates within the ",
      "scenario (1 = lowest; ties share the better rank). Rows are ordered by ",
      "true \\gls{pehe}. DR: the DR score. IF, true $\\tau$: the IF score ",
      "with the T-learner surrogate, which with the true nuisance functions ",
      "equals $\\tau$. IF, pseudo-outcome: the IF score with the \\gls{dr} ",
      "pseudo-outcome as its surrogate. MSE vs $\\mu_1$: the candidate's mean ",
      "squared distance from the true $\\mu_1$. Candidates: the true ",
      "\\gls{cate} $\\tau$; the \\gls{ate} $\\bar{\\tau}$; half-way between ",
      "them; $\\mu_1$; and $\\tau + X_2$. Table~\\ref{tab:me_if_selection} ",
      "makes the same comparison in the simulated trials."
    ),
    label = "me_if_population",
    landscape = FALSE,
    font_size = 10,
    tabcolsep = NULL
  )
  out_file <- file.path(tab_path, "me_if_population.tex")
  writeLines(tex, out_file)
  message("wrote ", out_file)
} else {
  message("skipped me_if_population.tex: no ", population_file,
          " (run `Rscript me_infl_investigation.R` in model_evaluation/)")
}

# ---- the IF score: selection by surrogate ------------------------------------

# per run, as me_selection() ranks: candidates with a finite score and truth,
# regret against the best of them. The correlations are across the candidates
# within the run; one is NA where a score ties across all of them.
if_runs <- metrics %>%
  pivot_longer(all_of(unname(IF_SEL_COLS)), names_to = "score",
               values_to = "value") %>%
  mutate(score = names(IF_SEL_COLS)[match(score, IF_SEL_COLS)]) %>%
  filter(is.finite(value), is.finite(true_pehe)) %>%
  group_by(scenario, n, run, score) %>%
  filter(n() >= 2) %>%
  summarise(
    regret = true_pehe[which.min(value)] - min(true_pehe),
    rho_pehe = suppressWarnings(cor(value, true_pehe, method = "spearman")),
    rho_mu1 = suppressWarnings(cor(value, mu1_mse, method = "spearman",
                                   use = "complete.obs")),
    .groups = "drop"
  )

# the DR and T-learner IF regrets are cells of me_regret.tex, from
# me_selection(); a disagreement means the two derivations have drifted
shared <- if_runs %>%
  filter(score %in% c("dr", "inflt")) %>%
  inner_join(
    selection %>%
      filter(paste0(score_type, "_", fold_type, "_", pipeline) %in%
               IF_SEL_COLS) %>%
      transmute(scenario, n, run, score = score_family, sel_regret = regret),
    by = c("scenario", "n", "run", "score")
  )
stopifnot(
  "me_if_selection's DR/IF regrets disagree with me_selection()'s" =
    nrow(shared) == sum(if_runs$score %in% c("dr", "inflt")) &&
    isTRUE(all.equal(shared$regret, shared$sel_regret))
)

if_wide <- if_runs %>%
  pivot_wider(id_cols = c(scenario, n, run), names_from = score,
              values_from = c(regret, rho_pehe, rho_mu1)) %>%
  full_join(
    ref_runs %>%
      filter(ref == "random") %>%
      select(scenario, n, run, regret_random = regret),
    by = c("scenario", "n", "run")
  )

if_summary <- summarise_metrics(
  if_wide, c("scenario", "n"),
  cols = c(dr = "regret_dr", infl = "regret_infl", inflt = "regret_inflt",
           random = "regret_random", rho_pehe = "rho_pehe_inflt",
           rho_mu1 = "rho_mu1_inflt"),
  count_na = character()
) %>%
  arrange(scenario, n)

if_regret_digits <- mcse_digits(unlist(
  if_summary[paste0("mcse_", c("dr", "infl", "inflt", "random"))]
))

if_cols <- tibble::tribble(
  ~stem,      ~header,                                       ~digits,
  "dr",       "DR",                                          if_regret_digits,
  "infl",     "\\shortstack[r]{IF,\\\\pseudo-outcome}",      if_regret_digits,
  "inflt",    "\\shortstack[r]{IF,\\\\T-learner}",           if_regret_digits,
  "random",   "\\shortstack[r]{Random\\\\selection}",        if_regret_digits,
  "rho_pehe", "\\shortstack[r]{True\\\\PEHE}",               2,
  "rho_mu1",  "\\shortstack[r]{MSE\\\\vs $\\mu_1$}",          2
)

if_cells <- lapply(seq_len(nrow(if_cols)), function(i) {
  fmt_mcse(if_summary[[paste0("mean_", if_cols$stem[i])]],
           if_summary[[paste0("mcse_", if_cols$stem[i])]],
           if_cols$digits[i])
})
names(if_cells) <- if_cols$stem

# a correlation missing from a run (a score tied across every candidate)
# leaves that run out of its mean only, so say how many
rho_na <- sum(is.na(if_wide$rho_pehe_inflt) | is.na(if_wide$rho_mu1_inflt))
rho_na_text <- if (rho_na > 0) {
  paste0(" The correlations leave out the ", rho_na, " simulated trial(s) ",
         "in which the score tied across every candidate.")
} else {
  ""
}

tex <- grouped_longtable(
  bind_cols(tibble::tibble(n = as.character(if_summary$n)),
            tibble::as_tibble(if_cells)),
  col_names = c("$n$", if_cols$header),
  group = if_summary$scenario,
  block = NULL,
  # short, or the spanning header widens the two columns under it; the
  # caption says whose correlation it is
  header_above = c(" " = 1, "Mean regret" = 4, "Spearman corr. with" = 2),
  align = c("l", rep("r", nrow(if_cols))),
  caption.short = "Model evaluation: the IF score's surrogate",
  caption = paste0(
    "The IF score with different surrogates in the simulated trials, with ",
    "the surrogate's nuisance functions fitted by XGBoost on the candidates' ",
    "training folds (CV) and the propensity score fixed at 0.5. Mean regret ",
    "of selecting by the DR score, by the IF score with the \\gls{dr} ",
    "pseudo-outcome as its surrogate, and by the IF score with the T-learner ",
    "surrogate $\\tilde{\\mu}_1 - \\tilde{\\mu}_0$, and of random selection. ",
    "Then the Spearman correlation, across the nine candidates within a ",
    "simulated trial, of the IF score (T-learner surrogate) with the ",
    "candidates' true \\gls{pehe} and with their mean squared distance from ",
    "the true $\\mu_1$. Mean (Monte Carlo SE) over ",
    runs_per_cell(if_runs, c("scenario", "n", "score")),
    " simulated trials per cell.", rho_na_text, " The DR and IF (T-learner) ",
    "columns repeat the CV, XGBoost, $\\pi$-fixed cells of ",
    "Table~\\ref{tab:me_regret}; Table~\\ref{tab:me_if_population} makes the ",
    "same comparison with the true nuisance functions."
  ),
  label = "me_if_selection",
  landscape = FALSE,
  font_size = 8,
  tabcolsep = "3pt"
)
out_file <- file.path(tab_path, "me_if_selection.tex")
writeLines(tex, out_file)
message("wrote ", out_file)
