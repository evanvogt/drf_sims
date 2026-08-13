##########
# title: figures for the MI confidence-interval example
##########
# The diagnostic counterpart to results_processing/thesis_figures/miss_cts_ci.R:
# that script picks the panels the chapter prints, this one shows every metric
# cts_miss_ci_metrics.R computes. Figures go to results/all_figures/, NOT
# results/thesis_figures/, so the two never overwrite each other.
#
# The factor structure differs from the other two missing-data studies:
# mechanism and method are single-valued (MAR, multiple_imputation), and what
# varies is scenario x model x STRATEGY - the three ways combine_mi_ci()
# (R/bootstrap_ci.R) pools the per-imputation bootstrap intervals. Four
# estimators only (CI_MODELS - no SuperLearner arm).

rm(list = ls())
# libraries
library(readr)
library(here)
library(patchwork)
library(purrr)
source(here("R", "figures.R"))

# paths
path <- here()
res_path <- file.path(dirname(path), "results", "missing", "ci_example")
fig_path <- file.path(dirname(path), "results", "all_figures", "missing", "ci_example")
dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)

metrics <- readRDS(file.path(res_path, "cts_miss_ci_metrics.RDS"))

# tidy up. Unlike the chapter script this keeps all five scenarios - this study's
# grid runs 1:5, and scenario 3 ("Two HTE vars") exists here even though
# missing/continuous and missing/binary never run it.
metrics <- metrics %>%
  apply_labels(MISS_SCENARIO_LABELS) %>%
  mutate(n = factor(n), prop = factor(prop)) %>%
  droplevels()

metrics_summary <- summarise_metrics(
  metrics,
  c("scenario", "n", "type", "prop", "mechanism", "method", "model", "strategy"),
  cols = c(bias = "bias", ate_bias = "ate_bias", mse = "mse", rmse = "rmse",
           mae = "mae", corr = "corr", spearman = "spearman",
           sign_acc = "sign_acc", mar_cov = "marginal_coverage",
           sim_cov = "simultaneous_coverage", ci_len = "mean_ci_length")
)

# helper: per-run distribution, strategy on x and one colour per estimator.
# Only one mechanism here, so facet on scenario alone.
ci_box_plot <- function(metrics, metric, y_lab, hline = NULL,
                        facet_scales = "free_y") {
  p <- ggplot(metrics, aes(x = strategy, y = .data[[metric]], colour = model)) +
    geom_boxplot(fill = "transparent", outlier.shape = NA) +
    facet_grid(cols = vars(scenario), scales = facet_scales) +
    drf_scale() +
    drf_theme(rotate_x = TRUE) +
    labs(y = y_lab, x = "Pooling strategy", colour = "Model")
  if (!is.null(hline)) {
    p <- p + geom_hline(yintercept = hline, linetype = "dashed")
  }
  p
}

# helper: the summary counterpart. point_range_plot()'s faceting is hardcoded to
# mechanism x scenario, which is degenerate here (one mechanism), so the facet is
# overridden at the call site - the same trick miss_cts_ci.R uses.
ci_sum_plot <- function(summary, metric, y_lab, hline = NULL,
                        facet_scales = "free_y") {
  p <- point_range_plot(summary, metric, y_lab, x = "strategy",
                        colour = "model", facet_scales = facet_scales) +
    facet_grid(cols = vars(scenario), scales = facet_scales) +
    labs(x = "Pooling strategy", colour = "Model")
  if (!is.null(hline)) {
    p <- p + geom_hline(yintercept = hline, linetype = "dashed")
  }
  p
}

# --- marginal coverage ------------------------------------------------------
# the proportion of units their own interval covers. Nominal 95%.
mar_cov_plot <- ci_box_plot(metrics, "marginal_coverage", "Marginal coverage",
                            hline = 0.95)
save_fig("cts_miss_ci_mar_cov_all.png", fig_path)

mar_cov_sum_plot <- ci_sum_plot(metrics_summary, "mar_cov", "Marginal coverage",
                                hline = 0.95)
save_fig("cts_miss_ci_mar_cov_summary.png", fig_path)

# --- simultaneous coverage --------------------------------------------------
# whether the band covers every unit at once - what the half-sample bootstrap is
# constructed to control. Evaluated per-unit on the sample's own covariates, so
# NOT comparable to the grid-based simultaneous coverage confidence_intervals/
# now reports (see this study's README).
sim_cov_plot <- ci_box_plot(metrics, "simultaneous_coverage",
                            "Simultaneous coverage", hline = 0.95)
save_fig("cts_miss_ci_sim_cov_all.png", fig_path)

sim_cov_sum_plot <- ci_sum_plot(metrics_summary, "sim_cov",
                                "Simultaneous coverage", hline = 0.95)
save_fig("cts_miss_ci_sim_cov_summary.png", fig_path)

# --- CI length --------------------------------------------------------------
# the width coverage is bought with. No reference line - there is no target
# width, only the coverage/width trade-off against the two panels above.
ci_len_plot <- ci_box_plot(metrics, "mean_ci_length", "CI length")
save_fig("cts_miss_ci_len_all.png", fig_path)

ci_len_sum_plot <- ci_sum_plot(metrics_summary, "ci_len", "CI length")
save_fig("cts_miss_ci_len_summary.png", fig_path)

# --- bias -------------------------------------------------------------------
# NOTE for every point-metric section below: the three strategies share one
# model_res$tau, so bias/MSE/RMSE/MAE/correlation/sign accuracy are IDENTICAL
# across strategies within a model. The triplicate points are the same number
# three times, not three disagreeing estimates - only the interval metrics above
# actually distinguish the pooling rules.
bias_plot <- ci_box_plot(metrics, "bias", "Bias", hline = 0)
save_fig("cts_miss_ci_bias_all.png", fig_path)

bias_sum_plot <- ci_sum_plot(metrics_summary, "bias", "Bias", hline = 0)
save_fig("cts_miss_ci_bias_summary.png", fig_path)

# --- ATE bias ---------------------------------------------------------------
ate_bias_plot <- ci_box_plot(metrics, "ate_bias", "ATE bias", hline = 0)
save_fig("cts_miss_ci_ate_bias_all.png", fig_path)

ate_bias_sum_plot <- ci_sum_plot(metrics_summary, "ate_bias", "ATE bias",
                                 hline = 0)
save_fig("cts_miss_ci_ate_bias_summary.png", fig_path)

# --- MSE --------------------------------------------------------------------
mse_plot <- ci_box_plot(metrics, "mse", "MSE", hline = 0)
save_fig("cts_miss_ci_mse_all.png", fig_path)

mse_sum_plot <- ci_sum_plot(metrics_summary, "mse", "MSE", hline = 0)
save_fig("cts_miss_ci_mse_summary.png", fig_path)

# --- RMSE -------------------------------------------------------------------
rmse_plot <- ci_box_plot(metrics, "rmse", "RMSE", hline = 0)
save_fig("cts_miss_ci_rmse_all.png", fig_path)

rmse_sum_plot <- ci_sum_plot(metrics_summary, "rmse", "RMSE", hline = 0)
save_fig("cts_miss_ci_rmse_summary.png", fig_path)

# --- MAE --------------------------------------------------------------------
mae_plot <- ci_box_plot(metrics, "mae", "MAE", hline = 0)
save_fig("cts_miss_ci_mae_all.png", fig_path)

mae_sum_plot <- ci_sum_plot(metrics_summary, "mae", "MAE", hline = 0)
save_fig("cts_miss_ci_mae_summary.png", fig_path)

# --- correlation with the truth ---------------------------------------------
# scenario 1 (Null) has no heterogeneity, so cate_metrics() reports corr and
# spearman as 0 rather than NA - excluded here. No reference line: 0 is not a
# target on this axis.
corr_metrics <- filter(metrics, scenario != "Null")
corr_summary <- filter(metrics_summary, scenario != "Null")

corr_plot <- ci_box_plot(corr_metrics, "corr", "Pearson correlation")
save_fig("cts_miss_ci_corr_all.png", fig_path)

corr_sum_plot <- ci_sum_plot(corr_summary, "corr", "Pearson correlation")
save_fig("cts_miss_ci_corr_summary.png", fig_path)

# --- Spearman rank correlation with the truth -------------------------------
spearman_plot <- ci_box_plot(corr_metrics, "spearman", "Spearman correlation")
save_fig("cts_miss_ci_spearman_all.png", fig_path)

spearman_sum_plot <- ci_sum_plot(corr_summary, "spearman",
                                 "Spearman correlation")
save_fig("cts_miss_ci_spearman_summary.png", fig_path)

# --- sign accuracy ----------------------------------------------------------
sign_acc_plot <- ci_box_plot(metrics, "sign_acc", "Proportion correct sign",
                             hline = 0.5)
save_fig("cts_miss_ci_sign_acc_all.png", fig_path)

sign_acc_sum_plot <- ci_sum_plot(metrics_summary, "sign_acc",
                                 "Proportion correct sign", hline = 0.5)
save_fig("cts_miss_ci_sign_acc_summary.png", fig_path)

# --- missing estimates (n_na) diagnostic ------------------------------------
# the count of NA CATE estimates within a run. A data-quality diagnostic, not a
# performance metric, and mostly zero - so a table listing only the cells that
# have any.
na_table <- metrics %>%
  group_by(scenario, model, strategy) %>%
  summarise(mean_n_na = mean(n_na, na.rm = TRUE),
            total_n_na = sum(n_na, na.rm = TRUE),
            .groups = "drop") %>%
  filter(total_n_na > 0) %>%
  arrange(scenario, model, strategy)

if (nrow(na_table) > 0) {
  print(na_table, n = Inf)
} else {
  print("no missing CATE estimates in any run")
}

# --- headline table ---------------------------------------------------------
headline <- metrics_summary %>%
  select(scenario, model, strategy, mean_mar_cov, mean_sim_cov, mean_ci_len,
         mean_bias, mean_mse, mean_rmse, mean_mae, mean_sign_acc, mean_corr) %>%
  arrange(scenario, model, strategy)

print(headline, n = Inf)
saveRDS(headline, file.path(res_path, "cts_miss_ci_headline.RDS"))
write_csv(headline, file.path(res_path, "cts_miss_ci_headline.csv"))

print("figures written!")
