##########
# title: LaTeX tables for the thesis appendix - the case-study learners
##########
# The supplementary model-evaluation analysis (model_evaluation/ADEMP.md, "The
# case-study learners"): the applied chapter's three learners - DR-RF, the
# causal forest (CF) and DR-SL - scored by the DR score in the model-evaluation
# simulated trials, with XGBoost as the surrogate nuisance learner on the CV
# folds. Separate from me_tables.R, whose tables it leaves untouched.
#
# - Per simulated trial, over the 3 learners: top-1 accuracy, regret, and for
#   each pair whether the DR score orders the two as their true PEHE does.
#   Regret's references are me_selection.R's, over the same 3.
# - Cells are "mean (Monte Carlo SE)" over simulated trials; 0/1 measures take
#   the binomial SE, as in me_tables.R.
# - The rank of DR-RF and CF among all 11 learners (the 9 candidates plus
#   these 2) is descriptive: the main study's selection is over 9.
#
# Writes to ../results/thesis_tables/: me_cs_learners, me_cs_ranks (.tex);
# with `extra`, also me_cs_learners_extra.tex, the estimated-pi version and the
# other two pairs. Without it they are printed to the console only, to judge
# whether they differ enough to report.
#
#   Rscript me_cs_tables.R [me_cs_metrics.RDS] [me_strat_metrics.RDS] [table dir] [extra]
#   Rscript me_cs_tables.R "" "" "" extra    # the defaults, plus the extra table

library(here)
source(here("R", "figures.R"))
source(here("R", "tables.R"))
source(here("model_evaluation", "me_selection.R"))

# paths
args <- commandArgs(trailingOnly = TRUE)
path <- here()
res <- file.path(dirname(path), "results")
# an absent or empty ("") argument takes the default
arg_or <- function(i, default) {
  if (length(args) >= i && nzchar(args[i])) args[i] else default
}
cs_file <- arg_or(1, file.path(res, "correlated", "model_evaluation_cs_learners",
                               "me_cs_metrics.RDS"))
strat_file <- arg_or(2, file.path(res, "correlated", "model_evaluation_strategies",
                                  "me_strat_metrics.RDS"))
tab_path <- arg_or(3, file.path(res, "thesis_tables"))
write_extra <- isTRUE(args[4] == "extra")
dir.create(tab_path, showWarnings = FALSE, recursive = TRUE)

# display labels, in column order
LEARNER_LABELS <- c(drrf = "DR-RF", cf = "CF", SL1 = "DR-SL")
PAIRS <- list(drrf_cf = c("drrf", "cf"), drrf_sl = c("drrf", "SL1"),
              cf_sl = c("cf", "SL1"))
PAIR_LABELS <- c(drrf_cf = "DR-RF vs CF", drrf_sl = "DR-RF vs DR-SL",
                 cf_sl = "CF vs DR-SL")
SCORE_LABELS <- c(dr05 = "$\\pi$-fixed", dr = "$\\pi$ estimated")

# ---- data --------------------------------------------------------------------

cs_raw <- readRDS(cs_file)
display <- function(d) {
  d %>%
    filter(scenario %in% names(SS_SCENARIO_LABELS)) %>%
    mutate(scenario = label_factor(scenario, SS_SCENARIO_LABELS),
           n = factor(n, levels = sort(unique(n))))
}
cs <- display(cs_raw)

# ---- selection among the 3, per simulated trial ------------------------------

#' 1 when the score orders learners a and b as their true PEHE does
pair_match <- function(model, score, pehe, a, b) {
  i <- match(c(a, b), model)
  as.numeric(sign(score[i[1]] - score[i[2]]) == sign(pehe[i[1]] - pehe[i[2]]))
}

#' One row per (scenario, n, run, score): hit, regret and the pair orderings
cs_selection <- function(cs, score) {
  cs %>%
    filter(is.finite(.data[[score]]), is.finite(true_pehe)) %>%
    group_by(scenario, n, run) %>%
    filter(n() == length(LEARNER_LABELS)) %>%
    summarise(
      hit = as.numeric(model[which.min(.data[[score]])] ==
                         model[which.min(true_pehe)]),
      regret = true_pehe[which.min(.data[[score]])] - min(true_pehe),
      !!!lapply(PAIRS, function(p) {
        rlang::expr(pair_match(model, .data[[!!score]], true_pehe, !!p[1], !!p[2]))
      }),
      .groups = "drop"
    ) %>%
    mutate(score = score)
}

selection <- bind_rows(lapply(names(SCORE_LABELS), cs_selection, cs = cs))

sel_summary <- summarise_metrics(
  selection, c("scenario", "n", "score"),
  cols = c(hit = "hit", regret = "regret",
           setNames(names(PAIRS), names(PAIRS))),
  binomial = c("hit", names(PAIRS)),
  count_na = character()
)

# regret references over the same 3 learners; they do not depend on the score
refs <- me_references(cs) %>%
  summarise_metrics(c("scenario", "n", "ref"), cols = c(regret = "regret"),
                    count_na = character()) %>%
  select(scenario, n, ref, mean_regret, mcse_regret) %>%
  tidyr::pivot_wider(names_from = ref, values_from = c(mean_regret, mcse_regret))

# mean true PEHE per learner, one column each
pehe_summary <- cs %>%
  select(scenario, n, run, model, true_pehe) %>%
  tidyr::pivot_wider(names_from = model, values_from = true_pehe,
                     names_prefix = "pehe_") %>%
  summarise_metrics(
    c("scenario", "n"),
    cols = setNames(paste0("pehe_", names(LEARNER_LABELS)),
                    paste0("pehe_", names(LEARNER_LABELS))),
    count_na = character()
  )

#' The summary columns of one score, renamed <stem>_<score>
score_wide <- function(score) {
  sel_summary %>%
    filter(.data$score == !!score) %>%
    select(-score) %>%
    rename_with(~ paste0(.x, "_", score), -c(scenario, n))
}

summary_all <- pehe_summary %>%
  left_join(score_wide("dr05"), by = c("scenario", "n")) %>%
  left_join(score_wide("dr"), by = c("scenario", "n")) %>%
  left_join(refs, by = c("scenario", "n")) %>%
  arrange(scenario, n)

trials <- runs_per_cell(distinct(cs, scenario, n, run), c("scenario", "n"))

# ---- the tables ----------------------------------------------------------------

#' A (scenario, n) x measure table: n as the only row label, under a bold
#' scenario header
#'
#' @param cols stem, header and digits per column, as ss_latex_table() takes
#' @param header_above spanning headers, n's column first
cs_table <- function(summary, cols, header_above, ...) {
  cells <- lapply(seq_len(nrow(cols)), function(i) {
    fmt_mcse(summary[[paste0("mean_", cols$stem[i])]],
             summary[[paste0("mcse_", cols$stem[i])]], cols$digits[i])
  })
  names(cells) <- cols$stem
  tab <- bind_cols(tibble(n = as.character(summary$n)), as_tibble(cells))
  grouped_longtable(
    tab,
    col_names = c("$n$", cols$header),
    group = summary$scenario,
    block = NULL,
    header_above = header_above,
    align = c("l", rep("r", nrow(cols))),
    ...
  )
}

write_tex <- function(tex, file) {
  out_file <- file.path(tab_path, file)
  writeLines(tex, out_file)
  message("wrote ", out_file)
}

PEHE_COLS <- tibble(
  stem = paste0("pehe_", names(LEARNER_LABELS)),
  header = unname(LEARNER_LABELS),
  digits = 3
)

PEHE_NOTE <- paste0(
  "the DR-learner with regression-forest nuisances and second stage (DR-RF), ",
  "the causal forest (CF) and the DR-learner with SuperLearner nuisances and ",
  "second stage (DR-SL, candidate SL1 of the main study), scored by the DR ",
  "score in the simulated trials of Tables~\\ref{tab:me_top1}--",
  "\\ref{tab:me_regret}, with XGBoost as the surrogate nuisance learner fit on ",
  "the CV folds. DR-RF and CF are the learners of the sample-size study, fit ",
  "to the same simulated trials. "
)

main_cols <- bind_rows(
  PEHE_COLS,
  tibble(stem = c("hit_dr05", "regret_dr05", "drrf_cf_dr05"),
         header = c("Top-1", "Regret", "DR-RF vs CF"),
         digits = c(2, 3, 2)),
  tibble(stem = c("regret_random", "regret_best_fixed"),
         header = c("Random", "Best fixed"),
         digits = 3)
)

tex <- cs_table(
  summary_all, main_cols,
  header_above = c(" " = 1, "Mean true PEHE" = 3,
                   "DR score, $\\\\pi$-fixed" = 3, "Regret references" = 2),
  caption.short = "Model evaluation: the case-study learners",
  caption = paste0(
    "Model evaluation, the case-study learners: ", PEHE_NOTE,
    "Mean true PEHE of each learner; then, selecting among the three by the ",
    "DR score with the propensity fixed at the true 0.5 ($\\pi$-fixed), ",
    "top-1 accuracy (the proportion of simulated trials in which the metric ",
    "selects the learner with the lowest true PEHE; random selection's is ",
    "1/3), regret (the true PEHE of the selected learner minus the lowest of ",
    "the three) and the DR-RF vs CF ordering accuracy (the proportion of ",
    "simulated trials in which the DR score orders DR-RF and CF as their true ",
    "PEHE does). Regret references: random selection, and the best fixed ",
    "choice, always deploying the learner with the lowest mean regret in that ",
    "scenario and $n$, chosen with hindsight. Mean (Monte Carlo SE) over ",
    trials, " simulated trials per cell; the SE of top-1 and ordering ",
    "accuracy is the binomial one."
  ),
  label = "me_cs_learners"
)
write_tex(tex, "me_cs_learners.tex")

# ---- the estimated-pi version and the other pairs -----------------------------

extra_cols <- tibble(
  stem = c("hit_dr", "regret_dr", "drrf_cf_dr", "drrf_sl_dr", "cf_sl_dr",
           "drrf_sl_dr05", "cf_sl_dr05"),
  header = c("Top-1", "Regret", unname(PAIR_LABELS), unname(PAIR_LABELS[-1])),
  digits = c(2, 3, 2, 2, 2, 2, 2)
)

# printed always, so whether they differ materially can be judged before
# deciding to write the table
cat("\nestimated pi vs pi-fixed, and the other two pairs (means):\n")
print(as.data.frame(summary_all %>%
  transmute(scenario, n,
            hit_dr05 = mean_hit_dr05, hit_dr = mean_hit_dr,
            regret_dr05 = mean_regret_dr05, regret_dr = mean_regret_dr,
            drrf_cf_dr05 = mean_drrf_cf_dr05, drrf_cf_dr = mean_drrf_cf_dr,
            drrf_sl_dr05 = mean_drrf_sl_dr05, drrf_sl_dr = mean_drrf_sl_dr,
            cf_sl_dr05 = mean_cf_sl_dr05, cf_sl_dr = mean_cf_sl_dr) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)

if (write_extra) {
  tex <- cs_table(
    summary_all, extra_cols,
    header_above = c(" " = 1, "DR score, $\\\\pi$ estimated" = 5,
                     "DR score, $\\\\pi$-fixed" = 2),
    caption.short = "Model evaluation: the case-study learners, further measures",
    caption = paste0(
      "Model evaluation, the case-study learners, further measures: the ",
      "learners and simulated trials of Table~\\ref{tab:me_cs_learners}. ",
      "Selecting among the three by the DR score with the scoring nuisance's ",
      "estimated propensity, trimmed to $[0.05, 0.95]$: top-1 accuracy, ",
      "regret and the ordering accuracy of each pair; and with the propensity ",
      "fixed at the true 0.5 ($\\pi$-fixed), the ordering accuracy of the two ",
      "pairs not in Table~\\ref{tab:me_cs_learners}. Mean (Monte Carlo SE) ",
      "over ", trials, " simulated trials per cell; the SE of top-1 and ",
      "ordering accuracy is the binomial one."
    ),
    label = "me_cs_learners_extra"
  )
  write_tex(tex, "me_cs_learners_extra.tex")
}

# ---- rank among all 11 ---------------------------------------------------------

# the 9 candidates' true PEHE from the strategies tree (the main tree's, which
# it carries unchanged), with DR-RF and CF added per simulated trial
nine <- readRDS(strat_file) %>%
  select(scenario, n, run, model, true_pehe)
eleven <- bind_rows(
  nine,
  cs_raw %>% filter(model %in% c("drrf", "cf")) %>%
    select(scenario, n, run, model, true_pehe)
) %>%
  filter(is.finite(true_pehe)) %>%
  group_by(scenario, n, run) %>%
  mutate(learners = n(), rank = rank(true_pehe, ties.method = "min")) %>%
  ungroup()

short <- eleven %>% distinct(scenario, n, run, learners) %>% filter(learners < 11)
if (nrow(short) > 0) {
  warning(nrow(short), " simulated trial(s) rank fewer than 11 learners",
          call. = FALSE)
}

rank_summary <- eleven %>%
  filter(model %in% c("drrf", "cf")) %>%
  select(scenario, n, run, model, rank) %>%
  tidyr::pivot_wider(names_from = model, values_from = rank,
                     names_prefix = "rank_") %>%
  display() %>%
  summarise_metrics(c("scenario", "n"),
                    cols = c(rank_drrf = "rank_drrf", rank_cf = "rank_cf"),
                    count_na = character()) %>%
  arrange(scenario, n)

tex <- cs_table(
  rank_summary,
  tibble(stem = c("rank_drrf", "rank_cf"), header = c("DR-RF", "CF"),
         digits = 2),
  header_above = c(" " = 1, "Mean true rank of 11" = 2),
  caption.short = "Model evaluation: the case-study learners among the candidates",
  caption = paste0(
    "Model evaluation, the case-study learners among the candidates: the ",
    "true rank of DR-RF and CF (Table~\\ref{tab:me_cs_learners}) by true PEHE ",
    "among the nine candidates of Table~\\ref{tab:me_fixed} and these two, ",
    "1 being the lowest, tied learners sharing the better rank. Descriptive ",
    "only: the metrics select among the nine. Mean (Monte Carlo SE) over ",
    runs_per_cell(distinct(eleven, scenario, n, run), c("scenario", "n")),
    " simulated trials per cell."
  ),
  label = "me_cs_ranks",
  landscape = FALSE
)
write_tex(tex, "me_cs_ranks.tex")
