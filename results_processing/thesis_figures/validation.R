##########
# title: figures for the thesis chapter - interim-analysis validation
##########
# The interim-analysis validation study (validation/continuous/), scenarios 1-4
# at rho = 0.5: three CATE estimators fit on the first `interim_prop` of a
# trial of 1000, and what they find there tested on the rest. The exploratory,
# one-scenario-per-render version of these figures is
# validation/continuous/cts_val_results.qmd. HC3 standard errors throughout.
#
# Every figure has stage 1's share of the trial on x and the scenarios across,
# and shows the mean over runs with a 95% Monte Carlo ribbon (binomial MCSE for
# a proportion of runs, sd / sqrt(runs) for a mean contrast):
#   - val_sg_direction.png: how often the top / bottom 10% responder subgroup
#     validates in stage 2, i.e. p < 0.05 with the W x subgroup estimate in the
#     direction stage 1 implied (the qmd's prop_sig_ref1). Dashed line in the
#     null: alpha / 2, the one-sided false-positive rate;
#   - val_sg_shrinkage.png: DR-RandomForest only, the subgroup contrast stage 1
#     claimed (ref1) against its true value in stage 2 (true2) - a winner's
#     curse is ref1 above true2;
#   - val_topvar.png: whether stage 1's most important covariate (x_top) is a
#     true modifier (top row) and whether it interacts with treatment in stage 2
#     by the spline-projection test (bottom row), per estimator and importance
#     measure. Dashed line in the null: the nominal 5%.
#
# Labels, palette, theme and figure sizing come from R/figures.R. This script
# carries only the paths and this study's filters.
#
# Writes to ../results/thesis_figures/validation/.
#
# Optional arguments, for a smoke run on mock metrics:
#   Rscript validation.R [metrics.RDS] [figure dir]

library(here)
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
  tibble(scenario = s,
         x_top = Filter(function(v) te_uses(sc_row, v), paste0("X", 1:5)))
})) %>%
  mutate(is_mod = TRUE)

# est / ref1 / true2 and the p-value of each subgroup test, one row per
# run x subgroup
DIR_PATTERN <- "^(p|est|ref1|true2)_(.+)$"
sg_long <- subgroups %>%
  rename(p_top = top_pval, p_bottom = bottom_pval) %>%
  select(scenario, model, interim_prop, run, matches(DIR_PATTERN)) %>%
  pivot_longer(matches(DIR_PATTERN), names_to = c(".value", "test"),
               names_pattern = DIR_PATTERN) %>%
  filter(!is.na(p)) %>%
  mutate(test = factor(c(top = "Top 10%", bottom = "Bottom 10%")[test],
                       levels = c("Top 10%", "Bottom 10%"))) %>%
  apply_labels(SS_SCENARIO_LABELS)

top_var <- top_var %>%
  left_join(modifiers, by = c("scenario", "x_top")) %>%
  mutate(is_mod = coalesce(is_mod, FALSE),
         measure = label_factor(measure, MEASURE_LABELS)) %>%
  apply_labels(SS_SCENARIO_LABELS)

#' A dashed reference line in the null scenario's panels only
#'
#' @param yintercept the line's height
#' @param levels the scenario factor's levels, so the line lands in its facet
#' @param ... further facet columns to pin it to, e.g. panel = "..."
null_ref <- function(yintercept, levels, ...) {
  tibble(scenario = factor(SS_SCENARIO_LABELS[["1"]], levels = levels),
         yintercept = yintercept, ...)
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
    aes(group = interaction(.data[[colour]], .data[[linetype]]),
        linetype = .data[[linetype]])
  }
  d <- mutate(d, lo = mean - Z * mcse, hi = mean + Z * mcse)
  if (proportion) d <- mutate(d, lo = pmax(lo, 0), hi = pmin(hi, 1))
  ggplot(d, aes(x = interim_prop, y = mean, colour = .data[[colour]],
                ymin = lo, ymax = hi)) +
    grp +
    geom_ribbon(aes(fill = .data[[colour]]), alpha = 0.15, colour = NA) +
    geom_line() +
    geom_point(size = 1.2) +
    # a tick at every interim point, numbered only at the ends and the middle
    scale_x_continuous(breaks = INTERIM_PROPS,
                       labels = ifelse(INTERIM_PROPS %in% X_LABELLED,
                                       as.character(INTERIM_PROPS),
                                       ""),
                       minor_breaks = NULL,
                       # room for 0.25 and 0.75, which otherwise run into the
                       # next panel's labels
                       expand = expansion(add = 0.04)) +
    drf_scale() +
    scale_fill_paletteer_d("rcartocolor::Safe") +
    drf_theme() +
    theme(legend.position = "bottom", legend.box = "vertical",
          panel.spacing.x = unit(0.8, "lines"))
}

# ---- subgroup validation rate ------------------------------------------------

sg_rate <- sg_long %>%
  group_by(scenario, test, model, interim_prop) %>%
  summarise(runs = n(),
            mean = mean(p < ALPHA & sign(est) == sign(ref1)),
            mcse = sqrt(mean * (1 - mean) / runs),
            .groups = "drop")

sg_rate_fig <- ribbon_lines(sg_rate, "model") +
  geom_hline(data = null_ref(ALPHA / 2, levels(sg_rate$scenario)),
             aes(yintercept = yintercept), linetype = "dashed") +
  facet_grid(rows = vars(test), cols = vars(scenario), scales = "free_y") +
  labs(x = X_LAB,
       y = "Proportion of runs validated\n(p < 0.05, stage 1's direction)",
       colour = "Model", fill = "Model")
save_fig("val_sg_direction.png", fig_path, height = FIG_HEIGHT,
         plot = sg_rate_fig)

# ---- subgroup contrast: stage 1's claim and the truth, DR-RandomForest ---------

sg_shrink <- sg_long %>%
  filter(model == MODEL_LABELS[["dr_random_forest"]]) %>%
  pivot_longer(c(ref1, true2), names_to = "quantity", values_to = "value") %>%
  mutate(quantity = label_factor(quantity,
                                 c(ref1 = "Stage 1 estimate (ref1)",
                                   true2 = "Stage 2 truth (true2)"))) %>%
  group_by(scenario, test, quantity, interim_prop) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            mcse = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
            .groups = "drop")

sg_shrink_fig <- ribbon_lines(sg_shrink, "quantity", proportion = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(rows = vars(test), cols = vars(scenario), scales = "free_y") +
  labs(x = X_LAB, y = "Mean contrast (subgroup − rest)",
       colour = NULL, fill = NULL)
save_fig("val_sg_shrinkage.png", fig_path, height = FIG_HEIGHT,
         plot = sg_shrink_fig)

# ---- top covariate: a true modifier, and does it replicate ---------------------

TOPVAR_PANELS <- c(modifier = "Nominates a true modifier",
                   proj = "Projection test p < 0.05")

topvar_rate <- bind_rows(
  top_var %>% transmute(scenario, model, measure, interim_prop,
                        panel = "modifier", hit = is_mod),
  top_var %>% filter(!is.na(p_proj)) %>%
    transmute(scenario, model, measure, interim_prop,
              panel = "proj", hit = p_proj < ALPHA)
) %>%
  group_by(scenario, panel, model, measure, interim_prop) %>%
  summarise(runs = n(),
            mean = mean(hit),
            mcse = sqrt(mean * (1 - mean) / runs),
            .groups = "drop") %>%
  mutate(panel = label_factor(panel, TOPVAR_PANELS))

topvar_fig <- ribbon_lines(topvar_rate, "model", linetype = "measure") +
  geom_hline(data = null_ref(ALPHA, levels(topvar_rate$scenario),
                             panel = factor(TOPVAR_PANELS[["proj"]],
                                            levels = TOPVAR_PANELS)),
             aes(yintercept = yintercept), linetype = "dashed") +
  facet_grid(rows = vars(panel), cols = vars(scenario), scales = "free_y") +
  # the ribbon would otherwise paint a grey box behind each linetype key
  guides(linetype = guide_legend(override.aes = list(fill = NA))) +
  labs(x = X_LAB, y = "Proportion of runs", colour = "Model", fill = "Model",
       linetype = "Importance measure")
save_fig("val_topvar.png", fig_path, height = FIG_HEIGHT, plot = topvar_fig)
