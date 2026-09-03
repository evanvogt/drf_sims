##########
# title: figures for the missing-covariate study, continuous outcome
##########
# The diagnostic counterpart to results_processing/thesis_figures/miss_cts.R:
# that script picks the two or three panels the chapter prints, this one shows
# every metric cts_miss_metrics.R computes. Figures go to results/all_figures/,
# NOT results/thesis_figures/, so the two never overwrite each other.
#
# Labels, palette, the mean +/- MCSE summary and the point-range panel all come
# from R/figures.R. This script carries the paths, the per-run boxplot shape and
# the section list.

rm(list = ls())
# libraries
library(readr)
library(here)
library(patchwork)
library(purrr)
source(here("R", "figures.R"))

# paths
path <- here()
res_path <- file.path(dirname(path), "results", "missing", "continuous")
fig_path <- file.path(dirname(path), "results", "all_figures", "missing", "continuous")
dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)

metrics <- readRDS(file.path(res_path, "cts_miss_metrics.RDS"))

# tidy up. n / type / prop are single-valued in this design (500 / both / 0.3),
# so nothing facets on them - but they stay factors, and stay in the grouping
# below, so the script survives a grid that gains a level.
metrics <- metrics %>%
  apply_labels(MISS_SCENARIO_LABELS) %>%
  mutate(n = factor(n), prop = factor(prop)) %>%
  droplevels()

# derived 0/1 reject-indicators for the power summaries below -
# summarise_metrics()'s binomial argument expects a real per-run column, not
# an inline expression
metrics <- metrics %>%
  mutate(BLP_reject = as.numeric(BLP_p < 0.05),
         indep_cate_reject = as.numeric(indep_cate < 0.05),
         indep_po_reject = as.numeric(indep_po < 0.05))

# per (scenario, mechanism, method, model) summaries. summarise_metrics() drops
# any entry of `cols` the tibble does not carry, so rel_efficiency being absent
# is not an error here.
metrics_summary <- summarise_metrics(
  metrics,
  c("scenario", "n", "type", "prop", "mechanism", "method", "model"),
  cols = c(bias = "bias", ate_bias = "ate_bias", mse = "mse", rmse = "rmse",
           mae = "mae", corr = "corr", spearman = "spearman",
           sign_acc = "sign_acc", BLP = "BLP_p", indep_cate = "indep_cate",
           indep_po = "indep_po", rel_eff = "rel_efficiency",
           rel_bias_complete = "rel_bias_complete",
           rel_ate_bias = "rel_ate_bias", rel_bias_cate = "rel_bias_cate",
           power_BLP = "BLP_reject", power_indep_cate = "indep_cate_reject",
           power_indep_po = "indep_po_reject"),
  binomial = c("power_BLP", "power_indep_cate", "power_indep_po")
)

# helper: the per-run distribution behind each summary panel. distribution_plot()
# in R/figures.R is hardwired to x = n / facet_wrap(~scenario), which is the
# sample-size studies' shape - this study has n fixed and two extra factors, so
# it needs its own. Faceting matches point_range_plot() so the boxplot and the
# summary panel stack cleanly. hline = NULL omits the reference line.
miss_box_plot <- function(metrics, metric, y_lab, hline = 0,
                          facet_scales = "free") {
  p <- ggplot(metrics, aes(x = model, y = .data[[metric]], colour = method)) +
    geom_boxplot(fill = "transparent", outlier.shape = NA) +
    facet_grid(rows = vars(mechanism), cols = vars(scenario),
               scales = facet_scales) +
    drf_scale() +
    drf_theme(rotate_x = TRUE) +
    labs(y = y_lab, x = "Model", colour = "Method")
  if (!is.null(hline)) {
    p <- p + geom_hline(yintercept = hline, linetype = "dashed")
  }
  p
}

# --- bias -------------------------------------------------------------------
bias_plot <- miss_box_plot(metrics, "bias", "Bias")
save_fig("cts_miss_bias_all.png", fig_path)

bias_sum_plot <- point_range_plot(metrics_summary, "bias", "Bias")
save_fig("cts_miss_bias_summary.png", fig_path)

# --- ATE bias ---------------------------------------------------------------
# bias of the average effect (mean(est) - mean(true)), as opposed to bias's
# per-unit average - the two can differ in sign or magnitude when errors don't
# average out evenly across units
ate_bias_plot <- miss_box_plot(metrics, "ate_bias", "ATE bias")
save_fig("cts_miss_ate_bias_all.png", fig_path)

ate_bias_sum_plot <- point_range_plot(metrics_summary, "ate_bias", "ATE bias")
save_fig("cts_miss_ate_bias_summary.png", fig_path)

# --- MSE --------------------------------------------------------------------
mse_plot <- miss_box_plot(metrics, "mse", "MSE", facet_scales = "free_y")
save_fig("cts_miss_mse_all.png", fig_path)

mse_sum_plot <- point_range_plot(metrics_summary, "mse", "MSE",
                                 facet_scales = "free_y")
save_fig("cts_miss_mse_summary.png", fig_path)

# --- RMSE -------------------------------------------------------------------
# sqrt(mse) - the same ranking as MSE, on the outcome's own scale
rmse_plot <- miss_box_plot(metrics, "rmse", "RMSE", facet_scales = "free_y")
save_fig("cts_miss_rmse_all.png", fig_path)

rmse_sum_plot <- point_range_plot(metrics_summary, "rmse", "RMSE",
                                  facet_scales = "free_y")
save_fig("cts_miss_rmse_summary.png", fig_path)

# --- MAE --------------------------------------------------------------------
# mean absolute error - a robustness check on MSE, less sensitive to outliers
mae_plot <- miss_box_plot(metrics, "mae", "MAE", facet_scales = "free_y")
save_fig("cts_miss_mae_all.png", fig_path)

mae_sum_plot <- point_range_plot(metrics_summary, "mae", "MAE",
                                 facet_scales = "free_y")
save_fig("cts_miss_mae_summary.png", fig_path)

# --- relative efficiency ----------------------------------------------------
# mse / mse of the same (scenario, mechanism, run, model) under complete_data.
# The reference arm sits at exactly 1 by construction; above 1 is the price of
# the missingness plus the handling method. This is the metric specific to the
# missing-data studies - nothing in binary/ or continuous/ has it.
if ("rel_efficiency" %in% names(metrics)) {
  rel_eff_plot <- miss_box_plot(metrics, "rel_efficiency",
                                "Relative efficiency (vs complete data)",
                                hline = 1, facet_scales = "free_y")
  save_fig("cts_miss_rel_eff_all.png", fig_path)

  rel_eff_sum_plot <- point_range_plot(metrics_summary, "rel_eff",
                                       "Relative efficiency (vs complete data)",
                                       facet_scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed")
  save_fig("cts_miss_rel_eff_summary.png", fig_path)
} else {
  warning("rel_efficiency absent from the metrics file - skipping those figures")
}

# --- relative bias (vs complete data) ----------------------------------------
# bias / bias of the same (scenario, mechanism, run, model) under complete_data.
# Unlike relative efficiency this is a SIGNED ratio: 1 means the same bias as
# complete data, a negative value means the handling method flipped the sign of
# the bias, and it blows up when bias_complete is near zero - read alongside the
# plain bias panel above, not on its own. This is an ARM-TO-ARM ratio (like
# relative efficiency), not a bias fraction - hence the hline = 1.
if ("rel_bias_complete" %in% names(metrics)) {
  rel_bias_complete_plot <- miss_box_plot(metrics, "rel_bias_complete",
                                          "Relative bias (vs complete data)",
                                          hline = 1, facet_scales = "free_y")
  save_fig("cts_miss_rel_bias_complete_all.png", fig_path)

  rel_bias_complete_sum_plot <- point_range_plot(metrics_summary, "rel_bias_complete",
                                                 "Relative bias (vs complete data)",
                                                 facet_scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed")
  save_fig("cts_miss_rel_bias_complete_summary.png", fig_path)
} else {
  warning("rel_bias_complete absent from the metrics file - skipping those figures")
}

# --- relative ATE bias (vs true) ---------------------------------------------
# (mean(est) - mean(true)) / mean(true): the ATE's bias as a fraction of the
# true ATE - the standard parameter-level "relative bias" from the simulation
# literature. This is a BIAS FRACTION, not an arm-to-arm ratio, so 0 (not 1) is
# the reference line, same convention as the plain bias/ATE bias panels.
if ("rel_ate_bias" %in% names(metrics)) {
  rel_ate_bias_plot <- miss_box_plot(metrics, "rel_ate_bias",
                                     "Relative ATE bias (vs true)",
                                     hline = 0, facet_scales = "free_y")
  save_fig("cts_miss_rel_ate_bias_all.png", fig_path)

  rel_ate_bias_sum_plot <- point_range_plot(metrics_summary, "rel_ate_bias",
                                            "Relative ATE bias (vs true)",
                                            facet_scales = "free_y") +
    geom_hline(yintercept = 0, linetype = "dashed")
  save_fig("cts_miss_rel_ate_bias_summary.png", fig_path)
} else {
  warning("rel_ate_bias absent from the metrics file - skipping those figures")
}

# --- relative bias of individual CATEs (vs true) ------------------------------
# mean((est - true) / true) per unit - the literal per-unit relative bias.
# Units where true == 0 are NA'd out upstream (R/metrics.R's cate_metrics())
# rather than contributing +/-Inf. Bias fraction, not a ratio - hline = 0.
if ("rel_bias_cate" %in% names(metrics)) {
  rel_bias_cate_plot <- miss_box_plot(metrics, "rel_bias_cate",
                                      "Relative bias (vs true CATE)",
                                      hline = 0, facet_scales = "free_y")
  save_fig("cts_miss_rel_bias_cate_all.png", fig_path)

  rel_bias_cate_sum_plot <- point_range_plot(metrics_summary, "rel_bias_cate",
                                             "Relative bias (vs true CATE)",
                                             facet_scales = "free_y") +
    geom_hline(yintercept = 0, linetype = "dashed")
  save_fig("cts_miss_rel_bias_cate_summary.png", fig_path)
} else {
  warning("rel_bias_cate absent from the metrics file - skipping those figures")
}

# --- correlation with the truth ---------------------------------------------
# scenario 1 (Null) has no heterogeneity, so cate_metrics() reports corr and
# spearman as 0 rather than NA - excluded here, the convention every results
# file in this repo uses. No reference line: 0 is not a target on this axis.
corr_metrics <- filter(metrics, scenario != "Null")
corr_summary <- filter(metrics_summary, scenario != "Null")

corr_plot <- miss_box_plot(corr_metrics, "corr", "Pearson correlation",
                           hline = NULL, facet_scales = "free_y")
save_fig("cts_miss_corr_all.png", fig_path)

corr_sum_plot <- point_range_plot(corr_summary, "corr", "Pearson correlation",
                                  facet_scales = "free_y")
save_fig("cts_miss_corr_summary.png", fig_path)

# --- Spearman rank correlation with the truth -------------------------------
# less sensitive to outlying CATE estimates and to monotone miscalibration than
# Pearson. Same Null exclusion.
spearman_plot <- miss_box_plot(corr_metrics, "spearman", "Spearman correlation",
                               hline = NULL, facet_scales = "free_y")
save_fig("cts_miss_spearman_all.png", fig_path)

spearman_sum_plot <- point_range_plot(corr_summary, "spearman",
                                      "Spearman correlation",
                                      facet_scales = "free_y")
save_fig("cts_miss_spearman_summary.png", fig_path)

# --- sign accuracy ----------------------------------------------------------
# proportion of units where sign(est) == sign(true) - direction detection, not
# magnitude. 0.5 is the chance level for a coin-flip on sign.
sign_acc_plot <- miss_box_plot(metrics, "sign_acc", "Proportion correct sign",
                               hline = 0.5, facet_scales = "free_x")
save_fig("cts_miss_sign_acc_all.png", fig_path)

sign_acc_sum_plot <- point_range_plot(metrics_summary, "sign_acc",
                                      "Proportion correct sign",
                                      facet_scales = "free_x") +
  geom_hline(yintercept = 0.5, linetype = "dashed")
save_fig("cts_miss_sign_acc_summary.png", fig_path)

# --- true CATE test evaluation -----------------------------------------
# BLP and independence tests run on the true CATE and true nuisances instead
# of an estimator's (see missing/continuous/README.md's "True-CATE HTE test
# evaluation") - the ceiling each model's own BLP/indep_cate test is chasing.
# No indep_po counterpart: dr_oracle's own indep_po already covers that. Bound
# in as an extra "model" x-axis category only where it's plotted (the BLP/
# indep_cate panels that follow, and the new power section below) - every
# other plot is unaffected. multiple_imputation rows are NA here, same as
# everywhere else in this study (see true_cate_test_row(), R/cate_models.R).
#
# BLP never touches X (run_blp_whole() takes Y, W, Y0.hat, W.hat, tau only),
# so it is identical across `method` for the same run - that flatness is
# expected, not a bug. indep_cate DOES depend on X, which genuinely differs by
# handling method, so its true-CATE series varies by method for real: it shows
# whether a method's reshaping of X alone (CATE held fixed at the truth) helps
# or hurts the independence test's power.
true_cate_tests <- readRDS(file.path(res_path, "cts_miss_true_cate_tests.RDS")) %>%
  apply_labels(MISS_SCENARIO_LABELS) %>%
  mutate(n = factor(n), prop = factor(prop)) %>%
  droplevels()

true_cate_summary <- true_cate_tests %>%
  group_by(scenario, n, type, prop, mechanism, method) %>%
  summarise(
    mean_BLP = mean(BLP_p, na.rm = TRUE),
    mcse_BLP = sd(BLP_p, na.rm = TRUE) / sqrt(sum(!is.na(BLP_p))),
    mean_indep_cate = mean(indep_cate, na.rm = TRUE),
    mcse_indep_cate = sd(indep_cate, na.rm = TRUE) / sqrt(sum(!is.na(indep_cate))),
    # named mean_power_*/mcse_power_*, not power_*/mcse_power_* - matching
    # summarise_metrics()'s own naming convention below, which always prefixes
    # mean_/mcse_ even for its binomial columns
    mean_power_BLP = mean(BLP_p < 0.05, na.rm = TRUE),
    mcse_power_BLP = sqrt(mean_power_BLP * (1 - mean_power_BLP) / sum(!is.na(BLP_p))),
    mean_power_indep_cate = mean(indep_cate < 0.05, na.rm = TRUE),
    mcse_power_indep_cate = sqrt(mean_power_indep_cate * (1 - mean_power_indep_cate) /
                                   sum(!is.na(indep_cate))),
    .groups = "drop"
  ) %>%
  mutate(model = "True CATE")

true_cate_raw <- true_cate_tests %>% mutate(model = "True CATE")

model_levels_true <- c(levels(metrics$model), "True CATE")

# row-bind the true-CATE reference series onto a per-run or summarised data
# frame, extending `model`'s levels so it becomes one more x-axis category
bind_true_cate <- function(df, true_df, extra_cols) {
  bind_rows(
    mutate(df, model = factor(as.character(model), levels = model_levels_true)),
    true_df %>%
      mutate(model = factor(model, levels = model_levels_true)) %>%
      select(scenario, n, type, prop, mechanism, method, model, all_of(extra_cols))
  )
}

# --- BLP test p-values ------------------------------------------------------
# dr_random_forest used to carry no BLP or independence test under
# profile = "missing" (PROFILES in R/cate_models.R), so its p-values were NA
# throughout. Fixed: the flag is now dr_rf_tests = TRUE, and the finished
# results were back-filled in place by R/patch_hte_tests.R rather than re-run.
# If this model is still NA here, the patch job has not been run over these
# results yet - check patch_status in check_all_studies.md.
BLP_plot <- miss_box_plot(bind_true_cate(metrics, true_cate_raw, "BLP_p"),
                          "BLP_p", "p-value", hline = 0.05, facet_scales = "free_x")
save_fig("cts_miss_blp_all.png", fig_path)

BLP_sum_plot <- point_range_plot(
  bind_true_cate(metrics_summary, true_cate_summary, c("mean_BLP", "mcse_BLP")),
  "BLP", "p-value", facet_scales = "free_x") +
  geom_hline(yintercept = 0.05, linetype = "dashed")
save_fig("cts_miss_blp_summary.png", fig_path)

# --- CATE permutation test p-values -----------------------------------------
indep_cate_plot <- miss_box_plot(bind_true_cate(metrics, true_cate_raw, "indep_cate"),
                                 "indep_cate", "p-value", hline = 0.05,
                                 facet_scales = "free_x")
save_fig("cts_miss_indep_cate_all.png", fig_path)

indep_cate_sum_plot <- point_range_plot(
  bind_true_cate(metrics_summary, true_cate_summary, c("mean_indep_cate", "mcse_indep_cate")),
  "indep_cate", "p-value", facet_scales = "free_x") +
  geom_hline(yintercept = 0.05, linetype = "dashed")
save_fig("cts_miss_indep_cate_summary.png", fig_path)

# --- PO permutation test p-values -------------------------------------------
indep_po_plot <- miss_box_plot(metrics, "indep_po", "p-value", hline = 0.05,
                               facet_scales = "free_x")
save_fig("cts_miss_indep_po_all.png", fig_path)

indep_po_sum_plot <- point_range_plot(metrics_summary, "indep_po", "p-value",
                                      facet_scales = "free_x") +
  geom_hline(yintercept = 0.05, linetype = "dashed")
save_fig("cts_miss_indep_po_summary.png", fig_path)

# --- HTE test power (proportion rejecting at alpha=0.05) --------------------
# new for this study - continuous/ and binary/ already report this; power_BLP/
# power_indep_cate/power_indep_po were added to summarise_metrics()'s cols
# above. The true-CATE reference series applies to BLP and indep_cate only,
# same as the p-value panels; indep_po has no true-CATE counterpart anywhere
# in this evaluation (dr_oracle's own indep_po already covers it).
BLP_power_plot <- point_range_plot(
  bind_true_cate(metrics_summary, true_cate_summary, c("mean_power_BLP", "mcse_power_BLP")),
  "power_BLP", "Power", facet_scales = "free_x") +
  geom_hline(yintercept = 0.05, linetype = "dashed")
save_fig("cts_miss_blp_power.png", fig_path)

indep_cate_power_plot <- point_range_plot(
  bind_true_cate(metrics_summary, true_cate_summary,
                c("mean_power_indep_cate", "mcse_power_indep_cate")),
  "power_indep_cate", "Power", facet_scales = "free_x") +
  geom_hline(yintercept = 0.05, linetype = "dashed")
save_fig("cts_miss_indep_cate_power.png", fig_path)

indep_po_power_plot <- point_range_plot(metrics_summary, "power_indep_po", "Power",
                                        facet_scales = "free_x") +
  geom_hline(yintercept = 0.05, linetype = "dashed")
save_fig("cts_miss_indep_po_power.png", fig_path)

# --- missing estimates (n_na) diagnostic ------------------------------------
# the count of NA CATE estimates within a run. A data-quality diagnostic, not a
# performance metric, and mostly zero - so a table rather than a boxplot, listing
# only the cells that have any. Printed, not saved, same as every other
# intermediate table here.
na_table <- metrics %>%
  group_by(scenario, mechanism, method, model) %>%
  summarise(mean_n_na = mean(n_na, na.rm = TRUE),
            total_n_na = sum(n_na, na.rm = TRUE),
            .groups = "drop") %>%
  filter(total_n_na > 0) %>%
  arrange(scenario, mechanism, method, model)

if (nrow(na_table) > 0) {
  print(na_table, n = Inf)
} else {
  print("no missing CATE estimates in any run")
}

# --- headline table ---------------------------------------------------------
# best model first within each (scenario, mechanism, method) cell
headline <- metrics_summary %>%
  select(scenario, mechanism, method, model, mean_bias, mean_mse, mean_rmse,
         mean_mae, mean_sign_acc, mean_corr, any_of("mean_rel_eff"),
         any_of("mean_rel_bias_complete"), any_of("mean_rel_ate_bias"),
         any_of("mean_rel_bias_cate"),
         mean_BLP, mean_indep_cate, mean_indep_po) %>%
  arrange(scenario, mechanism, method, mean_mse)

print(headline, n = Inf)
saveRDS(headline, file.path(res_path, "cts_miss_headline.RDS"))
write_csv(headline, file.path(res_path, "cts_miss_headline.csv"))

print("figures written!")
