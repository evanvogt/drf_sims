##########
# title: figures for the thesis chapter - interim-analysis validation
##########
# Thesis figures for the interim-analysis validation study (validation/continuous/,
# scenarios 1-4, rho = 0.5), stage 1's share of the trial on x. Means over runs
# with a 95% Monte Carlo ribbon; HC3 standard errors. The exploratory version is
# validation/val_results.qmd.
#
#   val_sg_direction  top / bottom 10% subgroup validates in stage 2 (p < 0.05,
#                     in stage 1's direction); null line alpha / 2
#   val_sg_shrinkage_<model>  each model's claimed contrast (ref1) against its
#                     true stage 2 value (true2): winner's curse is ref1 > true2
#   val_topvar        is stage 1's top covariate a true modifier, and does it
#                     interact in stage 2 (spline projection)? null line 5%.
#                     Grey dotted, bottom row: the one-stage comparator
#                     (cts_val_onestage.R, spline form, whole trial, Holm) -
#                     in the null, any covariate rejected; elsewhere, the
#                     smallest-p covariate is a true modifier and rejected
#
# The rates are unconditional: over every run, a run without a stage 2 test (no
# one predicted into the subgroup, or a leverage-1 row under HC3) counting as a
# miss. The null lines are therefore upper bounds. val_sg_shrinkage alone is
# over the tested runs, so its two contrasts average the same runs. <model> is
# the model's key in MODEL_LABELS, e.g. val_sg_shrinkage_dr_random_forest.
#
# Each is saved as <stem>.png (each panel its own y axis: compare by the labels,
# not the heights) and <stem>_shared_y.png, in ../results/thesis_figures/validation/.
#
#   Rscript validation.R [metrics.RDS] [figure dir] [onestage.RDS]   # mock smoke run

library(here)
library(ggh4x)
source(here("R", "figures.R"))
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
fig_path <- if (length(args) >= 2) {
  args[2]
} else {
  file.path(dirname(path), "results", "thesis_figures", "validation")
}
onestage_file <- if (length(args) >= 3) {
  args[3]
} else {
  file.path(study$res_path, "cts_val_onestage.RDS")
}
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

ALPHA <- HTE_ALPHA
Z <- qnorm(0.975)
MEASURE_LABELS <- c(tevim = "TE-VIM", shap = "TreeSHAP")
X_LAB <- "Proportion of participants in stage 1"
INTERIM_PROPS <- sort(unique(study$grid$interim_prop))
X_LABELLED <- c(0.25, 0.5, 0.75)
FIG_HEIGHT <- 14

# ---- data --------------------------------------------------------------------

metrics <- readRDS(metrics_file)

# HC3 rows only: x_top, est, ref1 and true2 repeat in the classical rows
subgroups <- metrics$subgroups %>% filter(se == "hc3")
top_var <- metrics$top_var %>% filter(se == "hc3")

# the true modifiers of each scenario, read off its te_expr as the qmd does
rho <- unique(study$grid$rho)
sc_table <- resolve_set(corr_set("continuous", rho))
modifiers <- bind_rows(lapply(sort(unique(study$grid$scenario)), function(s) {
  sc_row <- sc_table[sc_table$scenario == s, ]
  tibble(
    scenario = s,
    x_top = Filter(function(v) te_uses(sc_row, v), paste0("X", 1:5))
  )
})) %>%
  mutate(is_mod = TRUE)

# est / ref1 / true2 and the p-value of each subgroup test, one row per
# run x subgroup - every run, with p NA where there was no stage 2 test
DIR_PATTERN <- "^(p|est|ref1|true2)_(.+)$"
sg_long <- subgroups %>%
  rename(p_top = top_pval, p_bottom = bottom_pval) %>%
  select(scenario, model, interim_prop, run, matches(DIR_PATTERN)) %>%
  pivot_longer(
    matches(DIR_PATTERN),
    names_to = c(".value", "test"),
    names_pattern = DIR_PATTERN
  ) %>%
  mutate(
    test = factor(
      c(top = "Top 10%", bottom = "Bottom 10%")[test],
      levels = c("Top 10%", "Bottom 10%")
    )
  ) %>%
  apply_labels(SS_SCENARIO_LABELS)

top_var <- top_var %>%
  left_join(modifiers, by = c("scenario", "x_top")) %>%
  mutate(
    is_mod = coalesce(is_mod, FALSE),
    measure = label_factor(measure, MEASURE_LABELS)
  ) %>%
  apply_labels(SS_SCENARIO_LABELS)

#' A dashed reference line in the null scenario's panels only
#'
#' @param yintercept the line's height
#' @param levels the scenario factor's levels, so the line lands in its facet
#' @param ... further facet columns to pin it to, e.g. panel = "..."
null_ref <- function(yintercept, levels, ...) {
  tibble(
    scenario = factor(SS_SCENARIO_LABELS[["1"]], levels = levels),
    yintercept = yintercept,
    ...
  )
}

#' The scenarios across, `rows` down
#'
#' @param rows the column faceted down
#' @param independent TRUE for a y axis per panel, FALSE for one per row,
#'   shared across the scenarios
val_facets <- function(rows, independent = TRUE) {
  facet_grid2(
    rows = vars(.data[[rows]]),
    cols = vars(scenario),
    scales = "free_y",
    independent = if (independent) "y" else "none"
  )
}

#' Save an unfaceted figure twice: <stem>.png with a y axis per panel, and
#' <stem>_shared_y.png with one y axis per row across the scenarios
save_both <- function(fig, stem, rows) {
  save_fig(
    paste0(stem, ".png"),
    fig_path,
    height = FIG_HEIGHT,
    plot = fig + val_facets(rows)
  )
  save_fig(
    paste0(stem, "_shared_y.png"),
    fig_path,
    height = FIG_HEIGHT,
    plot = fig + val_facets(rows, independent = FALSE)
  )
}

#' Mean over runs against stage 1's share, a 95% Monte Carlo ribbon around it
#'
#' @param d summary with interim_prop, mean, mcse and the facet columns
#' @param colour column mapped to colour and to the ribbon's fill
#' @param linetype optional column mapped to linetype
#' @param proportion TRUE when `mean` is a proportion of runs: the ribbon is
#'   cut at 0 and 1, so it doesn't stretch the free y axis past either
ribbon_lines <- function(d, colour, linetype = NULL, proportion = TRUE) {
  grp <- if (is.null(linetype)) {
    aes(group = .data[[colour]])
  } else {
    aes(
      group = interaction(.data[[colour]], .data[[linetype]]),
      linetype = .data[[linetype]]
    )
  }
  d <- mutate(d, lo = mean - Z * mcse, hi = mean + Z * mcse)
  if (proportion) {
    d <- mutate(d, lo = pmax(lo, 0), hi = pmin(hi, 1))
  }
  # the null column's proportion panels span at least 0-0.1, so they show 0
  # under the reference line, and the all-zero "true modifier" panel gets a 0
  # to 0.1 axis rather than -0.05 to 0.05. Every other panel's axis follows its
  # data. No facet-row column here, so it lands in every row of the null column.
  y_floor <- if (proportion) {
    geom_blank(
      data = tibble(
        scenario = factor(
          SS_SCENARIO_LABELS[["1"]],
          levels = levels(d$scenario)
        ),
        interim_prop = min(d$interim_prop),
        mean = c(0, 0.1)
      ),
      aes(x = interim_prop, y = mean),
      inherit.aes = FALSE
    )
  }
  # any other proportion panel flat within 0.1 (e.g. every run at 1) gets a
  # 0.1-wide axis held against 0 or 1, not ggplot's +/- 0.05 around the value
  y_span <- if (proportion) {
    facet_cols <- setdiff(
      names(d),
      c("interim_prop", "mean", "mcse", "lo", "hi", "runs", colour, linetype)
    )
    flat <- d %>%
      group_by(across(all_of(facet_cols))) %>%
      summarise(lo = min(lo), hi = max(hi), .groups = "drop") %>%
      filter(hi - lo < 0.1) %>%
      mutate(lo = pmax(0, pmin(lo, hi - 0.1)), hi = pmax(hi, lo + 0.1)) %>%
      pivot_longer(c(lo, hi), values_to = "mean") %>%
      mutate(interim_prop = min(d$interim_prop))
    geom_blank(data = flat, aes(x = interim_prop, y = mean), inherit.aes = FALSE)
  }
  ggplot(
    d,
    aes(
      x = interim_prop,
      y = mean,
      colour = .data[[colour]],
      ymin = lo,
      ymax = hi
    )
  ) +
    grp +
    y_floor +
    y_span +
    geom_ribbon(aes(fill = .data[[colour]]), alpha = 0.15, colour = NA) +
    geom_line() +
    geom_point(size = 1.2) +
    # a tick at every interim point, numbered only at the ends and the middle
    scale_x_continuous(
      breaks = INTERIM_PROPS,
      labels = ifelse(
        INTERIM_PROPS %in% X_LABELLED,
        as.character(INTERIM_PROPS),
        ""
      ),
      minor_breaks = NULL,
      # room for 0.25 and 0.75, which otherwise run into the
      # next panel's labels
      expand = expansion(add = 0.04)
    ) +
    drf_scale() +
    scale_fill_paletteer_d("rcartocolor::Safe") +
    drf_theme() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      panel.spacing.x = unit(0.4, "lines")
    )
}

# ---- subgroup validation rate ------------------------------------------------

sg_rate <- sg_long %>%
  group_by(scenario, test, model, interim_prop) %>%
  summarise(
    runs = n(),
    # a run without a test is not validated, rather than dropped
    mean = mean(coalesce(p < ALPHA & sign(est) == sign(ref1), FALSE)),
    mcse = sqrt(mean * (1 - mean) / runs),
    .groups = "drop"
  )

sg_rate_fig <- ribbon_lines(sg_rate, "model") +
  geom_hline(
    data = null_ref(ALPHA / 2, levels(sg_rate$scenario)),
    aes(yintercept = yintercept),
    linetype = "dashed"
  ) +
  labs(
    x = X_LAB,
    y = "Proportion of runs validated\n(p < 0.05, stage 1's direction)",
    colour = "Model",
    fill = "Model"
  )
save_both(sg_rate_fig, "val_sg_direction", "test")

# ---- subgroup contrast: stage 1's claim and the truth, per model ---------------

# the tested runs only: ref1 exists without a stage 2 subgroup but true2 does
# not, so over every run the two lines would average different runs
sg_shrink <- sg_long %>%
  filter(!is.na(p)) %>%
  pivot_longer(c(ref1, true2), names_to = "quantity", values_to = "value") %>%
  mutate(
    quantity = label_factor(
      quantity,
      c(ref1 = "Stage 1 estimate", true2 = "Stage 2 truth")
    )
  ) %>%
  group_by(scenario, test, model, quantity, interim_prop) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    mcse = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
    .groups = "drop"
  )

# one figure per model, named by its key in MODEL_LABELS
for (m in levels(droplevels(sg_shrink$model))) {
  key <- coalesce(names(MODEL_LABELS)[match(m, MODEL_LABELS)], m)
  sg_shrink_fig <- sg_shrink %>%
    filter(model == m) %>%
    select(-model) %>%
    ribbon_lines("quantity", proportion = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      x = X_LAB,
      y = "Mean contrast (subgroup − rest)",
      colour = NULL,
      fill = NULL
    )
  save_both(sg_shrink_fig, paste0("val_sg_shrinkage_", key), "test")
}

# ---- top covariate: a true modifier, and does it replicate ---------------------

TOPVAR_PANELS <- c(
  modifier = "Nominates a true modifier",
  proj = "Projection test p < 0.05"
)

topvar_rate <- bind_rows(
  top_var %>%
    transmute(
      scenario,
      model,
      measure,
      interim_prop,
      panel = "modifier",
      hit = is_mod
    ),
  top_var %>%
    transmute(
      scenario,
      model,
      measure,
      interim_prop,
      panel = "proj",
      hit = coalesce(p_proj < ALPHA, FALSE)
    )
) %>%
  group_by(scenario, panel, model, measure, interim_prop) %>%
  summarise(
    runs = n(),
    mean = mean(hit),
    mcse = sqrt(mean * (1 - mean) / runs),
    .groups = "drop"
  ) %>%
  mutate(panel = label_factor(panel, TOPVAR_PANELS))

# the one-stage comparator, one value per scenario: the spline form, which the
# projection test matches. In the null any rejection (the family-wise error),
# elsewhere the smallest-p covariate a true modifier and rejected after Holm
onestage_ref <- if (file.exists(onestage_file)) {
  readRDS(onestage_file)$summary %>%
    filter(
      form == "spline",
      measure == if_else(scenario == 1, "any_rejected", "top_correct")
    ) %>%
    transmute(
      scenario = factor(
        SS_SCENARIO_LABELS[as.character(scenario)],
        levels = levels(topvar_rate$scenario)
      ),
      panel = factor(TOPVAR_PANELS[["proj"]], levels = TOPVAR_PANELS),
      yintercept = mean
    )
} else {
  message("no one-stage comparator at ", onestage_file,
          ": val_topvar drawn without its line")
  NULL
}

topvar_fig <- ribbon_lines(topvar_rate, "model", linetype = "measure") +
  geom_hline(
    data = null_ref(
      ALPHA,
      levels(topvar_rate$scenario),
      panel = factor(TOPVAR_PANELS[["proj"]], levels = TOPVAR_PANELS)
    ),
    aes(yintercept = yintercept),
    linetype = "dashed"
  ) +
  {
    if (!is.null(onestage_ref)) {
      geom_hline(
        data = onestage_ref,
        aes(yintercept = yintercept),
        linetype = "dotted",
        colour = "grey50"
      )
    }
  } +
  # the ribbon would otherwise paint a grey box behind each linetype key
  guides(linetype = guide_legend(override.aes = list(fill = NA))) +
  labs(
    x = X_LAB,
    y = "Proportion of runs",
    colour = "Model",
    fill = "Model",
    linetype = "Importance measure"
  )
save_both(topvar_fig, "val_topvar", "panel")
