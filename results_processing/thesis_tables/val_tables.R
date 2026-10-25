##########
# title: LaTeX tables for the thesis appendix - interim-analysis validation
##########
# Appendix tables for the interim-analysis validation study
# (validation/continuous/, scenarios 1-4, rho = 0.5), the counterpart of
# thesis_figures/validation.R.
#
# - A subgroup validates if p < 0.05 in stage 2 with the estimate in the
#   direction stage 1 implied (the sign of ref1).
# - HC3 standard errors everywhere but val_sg_se, which compares the two. Every
#   proportion is one 0/1 per run, so it takes the binomial MCSE.
# - Rates are unconditional: over every run in the cell, a run without a stage 2
#   test counting as not validated (not rejected). A test is missing when stage
#   1's tree predicts nobody into the subgroup in stage 2, or under HC3 when a
#   leverage-1 row loses the p-value; a message counts the second kind. The null
#   reference rates (alpha / 2, alpha) are therefore upper bounds.
# - In the nomination tables a covariate group is bold where it modifies the
#   effect in that scenario (te_uses()).
#
# Writes to ../results/thesis_tables/: val_sg_top, val_sg_bottom,
# val_sg_counts, val_nom_tevim, val_nom_shap, val_sg_se (.tex). \input{} with
# \usepackage{booktabs, longtable, pdflscape, array}.
#
#   Rscript val_tables.R [metrics.RDS] [table dir]   # mock smoke run

library(here)
source(here("R", "figures.R"))
source(here("R", "tables.R"))
source(here("validation", "continuous", "cts_val_config.R"))
source(here("validation", "continuous", "cts_val_dgms.R"))

# paths
args <- commandArgs(trailingOnly = TRUE)
path <- here()
metrics_file <- if (length(args) >= 1) {
  args[1]
} else {
  file.path(study$res_path, "cts_val_metrics.RDS")
}
tab_path <- if (length(args) >= 2) {
  args[2]
} else {
  file.path(dirname(path), "results", "thesis_tables")
}
dir.create(tab_path, showWarnings = FALSE, recursive = TRUE)

ALPHA <- HTE_ALPHA
DIGITS <- 2
LANDMARKS <- c(0.25, 0.5, 0.75)
# the null runs to 1000 at LANDMARKS, everything else 100 (cts_val_config.R)
NULL_RUNS <- max(study$grid$run[study$grid$scenario == 1])
MEASURE_LABELS <- c(tevim = "TE-VIM", shap = "TreeSHAP")
SE_LABELS <- c(classical = "Classical", hc3 = "HC3")
STAGE1_HEADER <- "Proportion of participants in stage 1"

# cell text is LaTeX; header_above text loses one backslash (R/tables.R)
SG_LABELS <- c(top = "Top 10\\%", bottom = "Bottom 10\\%")
SG_HEADERS <- c(top = "Top 10\\\\%", bottom = "Bottom 10\\\\%")

# x_top's groups, in column order, and the spanning header over each
NOM_GROUPS <- list(
  "X4" = "X4",
  "X3, X5" = c("X3", "X5"),
  "X1, X2" = c("X1", "X2"),
  "X01--X03" = c("X01", "X02", "X03"),
  "X04, X05" = c("X04", "X05")
)
NOM_HEADERS <- c("Modifier", "Other modifiers", "Prognostic",
                 "Correlated noise", "Independent noise")

interim_label <- function(x) formatC(x, format = "f", digits = 2)

#' "100", or "97--100" when the runs differ across the cells summarised
run_range <- function(runs) {
  r <- range(runs)
  if (r[1] == r[2]) as.character(r[1]) else paste0(r[1], "--", r[2])
}

write_table <- function(tex, name) {
  out_file <- file.path(tab_path, paste0(name, ".tex"))
  writeLines(tex, out_file)
  message("wrote ", out_file)
}

# ---- data --------------------------------------------------------------------

metrics <- readRDS(metrics_file)

# the true modifiers of each scenario, read off its te_expr as the figures do
rho <- unique(study$grid$rho)
sc_table <- resolve_set(corr_set("continuous", rho))
modifiers <- bind_rows(lapply(sort(unique(study$grid$scenario)), function(s) {
  sc_row <- sc_table[sc_table$scenario == s, ]
  tibble(scenario = s,
         x_top = Filter(function(v) te_uses(sc_row, v), paste0("X", 1:5)))
}))

# est / ref1 / true2 and the p-value of each subgroup test, one row per
# run x standard-error kind x subgroup. interim_prop is rounded so the
# landmarks match exactly
DIR_PATTERN <- "^(p|est|ref1|true2)_(.+)$"
sg_long <- metrics$subgroups %>%
  rename(p_top = top_pval, p_bottom = bottom_pval) %>%
  select(scenario, model, se, interim_prop, run, matches(DIR_PATTERN)) %>%
  pivot_longer(matches(DIR_PATTERN), names_to = c(".value", "test"),
               names_pattern = DIR_PATTERN) %>%
  mutate(interim_prop = round(interim_prop, 2),
         test = factor(test, levels = names(SG_LABELS))) %>%
  apply_labels(SS_SCENARIO_LABELS)

no_p <- sg_long %>% filter(!is.na(est), is.na(p)) %>% count(se)
if (nrow(no_p) > 0) {
  message("subgroup tests with a non-empty stage 2 subgroup but no p-value ",
          "(counted as not validated): ",
          paste0(no_p$se, " ", no_p$n, collapse = ", "))
}

# every run, tested or not: a missing p is NA, and coalesce() turns the NA a
# comparison with it gives into "not validated" rather than dropping the run
sg_hc3 <- filter(sg_long, se == "hc3")

# ---- subgroup validated in stage 1's direction, every interim point ----------

sg_rate <- sg_hc3 %>%
  mutate(validated = coalesce(p < ALPHA & sign(est) == sign(ref1), FALSE)) %>%
  summarise_metrics(c("test", "scenario", "model", "interim_prop"),
                    cols = c(validated = "validated"),
                    binomial = "validated", count_na = character()) %>%
  left_join(count(sg_hc3, test, scenario, model, interim_prop, name = "runs"),
            by = c("test", "scenario", "model", "interim_prop"))

sg_rate_caption <- function(s) {
  paste0(
    "Interim-analysis validation, ", SG_LABELS[[s]], " responder subgroup: ",
    "proportion of runs in which the subgroup validates in stage 2, mean ",
    "(Monte Carlo SE). Stage 1's estimated CATEs define its top and bottom ",
    "10\\% of responders; a classification tree learns them from the ",
    "covariates and predicts the subgroup into stage 2, where the ",
    "treatment-by-subgroup interaction is tested (HC3 standard errors). A run ",
    "validates if $p < ", ALPHA, "$ with the interaction in the direction ",
    "stage 1 implied; a run in which no stage 2 participant was predicted ",
    "into the subgroup has no test and counts as not validated. Runs: 100 (",
    NULL_RUNS, " at 0.25, 0.50 and 0.75 in the null scenario), as the range ",
    "over the interim points. In the null scenario every validation is a ",
    "false positive, at a rate of at most ", ALPHA / 2, "."
  )
}

for (s in names(SG_LABELS)) {
  d <- filter(sg_rate, test == s)
  interim_cols <- interim_label(sort(unique(d$interim_prop)))

  body <- d %>%
    transmute(scenario, model, interim = interim_label(interim_prop),
              cell = fmt_mcse(mean_validated, mcse_validated, DIGITS)) %>%
    pivot_wider(names_from = interim, values_from = cell, names_sort = TRUE) %>%
    left_join(d %>% group_by(scenario, model) %>%
                summarise(runs = run_range(runs), .groups = "drop"),
              by = c("scenario", "model")) %>%
    arrange(scenario, model)

  tab <- body %>%
    transmute(model = as.character(model), runs,
              across(all_of(interim_cols), ~ coalesce(.x, "---")))

  name <- paste0("val_sg_", s)
  write_table(grouped_longtable(
    tab,
    col_names = c("Estimator", "Runs", interim_cols),
    group = body$scenario,
    block = NULL,
    header_above = c(" " = 2, setNames(length(interim_cols), STAGE1_HEADER)),
    align = c("l", "r", rep("r", length(interim_cols))),
    caption = sg_rate_caption(s),
    caption.short = paste0("Interim-analysis validation: ", SG_LABELS[[s]],
                           " subgroup validated"),
    label = name
  ), name)
}

# ---- direction of the subgroup tests, landmark interim points ----------------

# true2 is 0 in the null (within rounding), where it has no direction
COUNT_STATS <- c("runs", "sig_same", "sig_opp", "true_same")

sg_counts <- sg_hc3 %>%
  filter(interim_prop %in% LANDMARKS) %>%
  mutate(sig = coalesce(p < ALPHA, FALSE),
         same = sign(est) == sign(ref1),
         true_same = ifelse(abs(true2) < 1e-8, NA, sign(true2) == sign(ref1))) %>%
  group_by(scenario, model, interim_prop, test) %>%
  summarise(runs = n(),
            sig_same = sum(sig & same),
            sig_opp = sum(sig & !same),
            true_same = if (all(is.na(true_same))) NA_integer_ else
              sum(true_same, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(across(all_of(COUNT_STATS), ~ ifelse(is.na(.x), "---", as.character(.x))))

count_cols <- paste(rep(names(SG_LABELS), each = length(COUNT_STATS)),
                    COUNT_STATS, sep = "__")
body <- sg_counts %>%
  pivot_longer(all_of(COUNT_STATS), names_to = "stat", values_to = "cell") %>%
  mutate(key = paste(test, stat, sep = "__")) %>%
  select(-test, -stat) %>%
  pivot_wider(names_from = key, values_from = cell) %>%
  arrange(scenario, model, interim_prop)
for (cc in setdiff(count_cols, names(body))) body[[cc]] <- NA_character_

tab <- body %>%
  transmute(model = as.character(model), interim = interim_label(interim_prop),
            across(all_of(count_cols), ~ coalesce(.x, "---")))

write_table(grouped_longtable(
  tab,
  col_names = c("Estimator", "Stage 1",
                rep(c("Runs", "Same", "Opposite", "True same"), length(SG_LABELS))),
  group = body$scenario,
  block = body$model,
  header_above = c(" " = 2, setNames(rep(length(COUNT_STATS), length(SG_HEADERS)),
                                     SG_HEADERS)),
  caption = paste0(
    "Interim-analysis validation: direction of the stage 2 subgroup tests, ",
    "number of runs, at three interim points (stage 1: the proportion of ",
    "participants in stage 1). Runs: every run, 100 (", NULL_RUNS, " in the ",
    "null scenario), including those in which no stage 2 participant was ",
    "predicted into the subgroup and so had no test, which count in neither ",
    "Same nor Opposite. Same: $p < ", ALPHA, "$ (HC3) with the treatment-by-subgroup ",
    "interaction in the direction stage 1's estimated CATEs implied (positive ",
    "for the top 10\\%, negative for the bottom 10\\%). Opposite: $p < ",
    ALPHA, "$ in the other direction. True same: the same contrast in stage ",
    "2's true CATEs has the direction stage 1 implied, whether or not the ",
    "test is significant; it has no direction in the null scenario, where ",
    "the true CATE is constant. In the null, Same and Opposite are each ",
    "expected at no more than ", ALPHA / 2, " of Runs."
  ),
  caption.short = "Interim-analysis validation: direction of the subgroup tests",
  label = "val_sg_counts",
  landscape = FALSE
), "val_sg_counts")

# ---- the covariate stage 1 nominates, landmark interim points ----------------

# x_top repeats in the classical rows
nom_group <- setNames(rep(names(NOM_GROUPS), lengths(NOM_GROUPS)),
                      unlist(NOM_GROUPS))
top_var <- metrics$top_var %>%
  filter(se == "hc3", round(interim_prop, 2) %in% LANDMARKS) %>%
  mutate(interim_prop = round(interim_prop, 2))
unknown <- setdiff(unique(top_var$x_top), names(nom_group))
if (length(unknown) > 0) {
  stop("x_top outside the covariate groups: ", paste(unknown, collapse = ", "))
}

# one 0/1 per run and group: did stage 1 nominate a covariate in it
nominations <- top_var %>%
  select(scenario, model, measure, interim_prop, run, x_top) %>%
  tidyr::crossing(group = names(NOM_GROUPS)) %>%
  mutate(nominated = nom_group[x_top] == group) %>%
  summarise_metrics(c("measure", "scenario", "model", "interim_prop", "group"),
                    cols = c(nom = "nominated"), binomial = "nom",
                    count_na = character())

# bold where the group modifies the effect in that scenario
mod_groups <- modifiers %>%
  transmute(scenario, group = nom_group[x_top], is_mod = TRUE) %>%
  distinct()
nominations <- nominations %>%
  left_join(mod_groups, by = c("scenario", "group")) %>%
  mutate(cell = fmt_mcse(mean_nom, mcse_nom, DIGITS),
         cell = ifelse(coalesce(is_mod, FALSE), paste0("\\textbf{", cell, "}"), cell)) %>%
  apply_labels(SS_SCENARIO_LABELS)

nom_counts <- count(top_var, measure, scenario, model, interim_prop)
nom_runs <- run_range(nom_counts$n[nom_counts$scenario != 1])
nom_runs_null <- run_range(nom_counts$n[nom_counts$scenario == 1])

for (m in names(MEASURE_LABELS)) {
  body <- nominations %>%
    filter(measure == m) %>%
    select(scenario, model, interim_prop, group, cell) %>%
    pivot_wider(names_from = group, values_from = cell) %>%
    arrange(scenario, model, interim_prop)

  tab <- body %>%
    transmute(model = as.character(model), interim = interim_label(interim_prop),
              across(all_of(names(NOM_GROUPS)), ~ coalesce(.x, "---")))

  name <- paste0("val_nom_", m)
  write_table(grouped_longtable(
    tab,
    col_names = c("Estimator", "Stage 1", names(NOM_GROUPS)),
    group = body$scenario,
    block = body$model,
    header_above = c(" " = 2, setNames(rep(1, length(NOM_HEADERS)), NOM_HEADERS)),
    caption = paste0(
      "Interim-analysis validation: the covariate stage 1 ranks most ",
      "important by ", MEASURE_LABELS[[m]], ", proportion of runs nominating ",
      "one in each group, mean (Monte Carlo SE) over ", nom_runs, " runs (",
      nom_runs_null, " in the null scenario), at ",
      "three interim points (stage 1: the proportion of participants in stage ",
      "1). Each row sums to 1. X4 modifies the effect in every scenario but ",
      "the null, X3 and X5 in the complex scenario only; elsewhere they are ",
      "noise correlated with X4. Bold: a group that modifies the effect in ",
      "that scenario. X1 and X2 are prognostic, X01--X03 noise correlated ",
      "with the other covariates, and X04 and X05, the indicators of one ",
      "three-level covariate, independent noise."
    ),
    caption.short = paste0("Interim-analysis validation: covariate nominated by ",
                           MEASURE_LABELS[[m]]),
    label = name,
    landscape = FALSE
  ), name)
}

# ---- subgroup tests, classical vs HC3 standard errors ------------------------

# over every run, so both kinds share a denominator: a leverage-1 run has a
# classical p but no HC3 one, and counts as not rejected under HC3
se_rate <- sg_long %>%
  filter(interim_prop %in% LANDMARKS) %>%
  mutate(rejected = coalesce(p < ALPHA, FALSE)) %>%
  summarise_metrics(c("scenario", "test", "model", "interim_prop", "se"),
                    cols = c(rej = "rejected"), binomial = "rej",
                    count_na = character())

se_cols <- paste(rep(interim_label(LANDMARKS), each = length(SE_LABELS)),
                 names(SE_LABELS), sep = "__")
body <- se_rate %>%
  transmute(scenario, test, model,
            key = paste(interim_label(interim_prop), se, sep = "__"),
            cell = fmt_mcse(mean_rej, mcse_rej, DIGITS)) %>%
  pivot_wider(names_from = key, values_from = cell) %>%
  arrange(scenario, test, model)
for (cc in setdiff(se_cols, names(body))) body[[cc]] <- NA_character_

tab <- body %>%
  transmute(test = unname(SG_LABELS[as.character(test)]),
            model = as.character(model),
            across(all_of(se_cols), ~ coalesce(.x, "---")))

write_table(grouped_longtable(
  tab,
  col_names = c("Subgroup", "Estimator",
                rep(unname(SE_LABELS), length(LANDMARKS))),
  group = body$scenario,
  block = body$test,
  header_above = list(
    c(" " = 2, setNames(rep(length(SE_LABELS), length(LANDMARKS)),
                        interim_label(LANDMARKS))),
    c(" " = 2, setNames(length(se_cols), STAGE1_HEADER))
  ),
  caption = paste0(
    "Interim-analysis validation: the stage 2 subgroup tests under classical ",
    "and HC3 standard errors, proportion of runs with $p < ", ALPHA, "$ in ",
    "either direction, mean (Monte Carlo SE), at three interim points. Over ",
    "every run; a run with no test (no stage 2 participant predicted into the ",
    "subgroup, or under HC3 a leverage-one observation) counts as not ",
    "rejected. In the null scenario this is the type I error rate, at most ",
    ALPHA, "."
  ),
  caption.short = "Interim-analysis validation: classical vs HC3 standard errors",
  label = "val_sg_se",
  landscape = FALSE
), "val_sg_se")
