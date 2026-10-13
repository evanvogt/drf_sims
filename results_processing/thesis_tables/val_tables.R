##########
# title: LaTeX tables for the thesis appendix - interim-analysis validation
##########
# The interim-analysis validation study (validation/continuous/), scenarios 1-4
# at rho = 0.5, the appendix counterpart of thesis_figures/validation.R: three
# CATE estimators fit on the first `interim_prop` of a trial of 1000, and what
# they find there tested on the rest.
#
# Tables:
#   - val_sg_top.tex, val_sg_bottom.tex: how often the top / bottom 10%
#     responder subgroup validates in stage 2 - p < 0.05 with the W x subgroup
#     estimate in the direction stage 1 implied (the sign of ref1), the
#     figure's val_sg_direction - at all 11 interim points. Rows scenario /
#     estimator, one table per subgroup so the cells keep their MCSE.
#     Landscape;
#   - val_sg_counts.tex: at the landmark interim points, per subgroup, the
#     number of runs significant in stage 1's direction, significant the other
#     way, and whose true stage 2 contrast (true2) has the sign stage 1 implied
#     - whether stage 1's claim was right, significant or not. Counts out of
#     Runs. Rows scenario / estimator / interim point;
#   - val_nom_tevim.tex, val_nom_shap.tex: which covariate stage 1 ranks most
#     important (x_top, the one the top-covariate tests are run on), grouped:
#     X4 | X3, X5 | X1, X2 (prognostic) | X01-X03 (noise correlated with the
#     modifiers) | X04, X05 (independent noise, the dummies of one three-level
#     categorical). At the landmarks, rows scenario / estimator / interim
#     point, one table per importance measure. A group that modifies the
#     effect in that scenario is bold: X4 in scenarios 2-4, X3 and X5 in
#     scenario 3 only (te_uses(), as the figures read it);
#   - val_sg_se.tex: the subgroup tests' rejection rate (p < 0.05, either
#     direction) under classical and HC3 standard errors, at the landmarks.
#     Rows scenario / subgroup / estimator.
#
# HC3 standard errors everywhere but val_sg_se. Every proportion is over runs,
# one 0/1 per run, so it gets the binomial MCSE, as the figures' ribbons do.
#
# Runs: a subgroup test exists only where stage 1's tree predicts at least one
# stage 2 participant into the subgroup (interaction_test(),
# validation/val_common.R); the runs without one are dropped from that cell,
# as in the figure. Under HC3 a test can also lack a p-value (a leverage-1
# row) though its subgroup is not empty - a message says how many runs that
# drops besides, since the captions only name the empty subgroups.
#
# Writes to ../results/thesis_tables/:
#   val_sg_top.tex, val_sg_bottom.tex, val_sg_counts.tex,
#   val_nom_tevim.tex, val_nom_shap.tex, val_sg_se.tex
# \input{} them into a document with
#   \usepackage{booktabs, longtable, pdflscape, array}
#
# Optional arguments, for a smoke run on mock metrics:
#   Rscript val_tables.R [metrics.RDS] [table dir]

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
          "(dropped as well): ",
          paste0(no_p$se, " ", no_p$n, collapse = ", "))
}

# the runs with a stage 2 test
sg_tested <- filter(sg_long, !is.na(p))
sg_hc3 <- filter(sg_tested, se == "hc3")

# ---- subgroup validated in stage 1's direction, every interim point ----------

sg_rate <- sg_hc3 %>%
  mutate(validated = p < ALPHA & sign(est) == sign(ref1)) %>%
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
    "stage 1 implied. Runs: the runs left after dropping those in which no ",
    "stage 2 participant was predicted into the subgroup, out of 100, as the ",
    "range over the interim points. In the null scenario every validation is ",
    "a false positive, at an expected rate of ", ALPHA / 2, "."
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
  mutate(sig = p < ALPHA,
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
    "participants in stage 1). Runs: the runs left after dropping those in ",
    "which no stage 2 participant was predicted into the subgroup, out of ",
    "100. Same: $p < ", ALPHA, "$ (HC3) with the treatment-by-subgroup ",
    "interaction in the direction stage 1's estimated CATEs implied (positive ",
    "for the top 10\\%, negative for the bottom 10\\%). Opposite: $p < ",
    ALPHA, "$ in the other direction. True same: the same contrast in stage ",
    "2's true CATEs has the direction stage 1 implied, whether or not the ",
    "test is significant; it has no direction in the null scenario, where ",
    "the true CATE is constant. In the null, Same and Opposite are each ",
    "expected at ", ALPHA / 2, " of Runs."
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

nom_runs <- run_range(count(top_var, measure, scenario, model, interim_prop)$n)

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
      "one in each group, mean (Monte Carlo SE) over ", nom_runs, " runs, at ",
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

se_rate <- sg_tested %>%
  filter(interim_prop %in% LANDMARKS) %>%
  mutate(rejected = p < ALPHA) %>%
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
    "the runs in which at least one stage 2 participant was predicted into ",
    "the subgroup. In the null scenario this is the type I error rate, ",
    "nominally ", ALPHA, "."
  ),
  caption.short = "Interim-analysis validation: classical vs HC3 standard errors",
  label = "val_sg_se",
  landscape = FALSE
), "val_sg_se")
