##########
# title: metrics for missing covariates, binary outcome
##########
# The pipeline lives in R/metrics.R. Grouping columns come from the study
# config, so this script carries only what is specific to missing covariates, binary outcome.

library(here)
source(here("missing/binary/bin_miss_config.R"))
source(here("R", "metrics.R"))
source(here("R", "cate_models.R"))
source(here("R", "missingness.R"))

# No bug N repair here any more: runs from the risk-difference DGM save the
# MNAR-Y (now MNAR-tau) truth already averaged over U (bU * tanh(U) has mean zero), and the
# logit-scale runs it repaired are superseded.
all_results_df <- readRDS(file.path(study$res_path, "bin_miss_all.RDS"))

# ... so none may be in the collection. Every logit-scale run made after the
# bug N fix carries truth$tau_u0, which no risk-difference run does.
has_tau_u0 <- vapply(all_results_df$results, function(rs) {
  any(vapply(rs, function(r) "tau_u0" %in% names(r$result$truth), logical(1)))
}, logical(1))
if (any(has_tau_u0)) {
  stop("bin_miss_all.RDS holds logit-scale results (truth$tau_u0) in ",
       sum(has_tau_u0), " parameter combination(s). Archive the old results ",
       "(R/archive_old_results.R) and re-collect.", call. = FALSE)
}

# MNAR-tau's secondary truth, tau + E[U_term | complete or incomplete] - see
# cate_metrics_split() in R/metrics.R and missing/ADEMP.md, "Estimands"
U_SHIFT <- u_term_by_completeness("binary_missing")

# every point metric on all analysed units, the complete units (_cu, the subset
# every method shares) and the incomplete units (_iu) - see cate_metrics_split()
metrics <- compute_metrics(
  study, all_results_df, models = CATE_MODELS,
  per_model = function(model_res, true_tau, model, sim_res, keys) {
    bind_cols(
      cate_metrics_split(model_res$tau, true_tau, keys$scenario,
                         sim_res$miss_mask, sim_res$retained_indices,
                         r_shift = if (keys$mechanism == "MNAR-tau") U_SHIFT),
      hte_test_metrics(model_res, sim_res, model)
    )
  }
)

# Comparisons against the complete-data reference arm, same (scenario,
# mechanism, run, model). The headline is rel_efficiency_cu, on the complete
# units, which every method has. The all-unit rel_efficiency is kept only for
# the methods that analyse all 500 units: for complete_cases and IPW it would
# divide a ~350-unit MSE by a 500-unit one, so it is NA there. The bias
# comparisons are differences, not ratios - bias_complete is often near zero,
# which made the old per-run ratio rel_bias_complete unstable (see
# missing/README.md). rel_ate_bias and rel_bias_cate (relative to the TRUE
# parameter) arrive from cate_metrics() in R/metrics.R; rel_bias_cate is not
# plotted here, since the true CATE crosses zero.
ref_by <- setdiff(c(study$path_cols, "run", "model"), "method")
ROW_DROPPING <- c("complete_cases", "IPW")

complete_ref <- metrics %>%
  filter(method == "complete_data") %>%
  select(all_of(ref_by), mse_complete = mse, mse_cu_complete = mse_cu,
         mse_iu_complete = mse_iu, bias_complete = bias,
         bias_cu_complete = bias_cu)

metrics <- metrics %>%
  left_join(complete_ref, by = ref_by) %>%
  mutate(rel_efficiency = if_else(method %in% ROW_DROPPING, NA_real_,
                                  mse / mse_complete),
         rel_efficiency_cu = mse_cu / mse_cu_complete,
         rel_efficiency_iu = mse_iu / mse_iu_complete,
         bias_diff_complete = if_else(method %in% ROW_DROPPING, NA_real_,
                                      bias - bias_complete),
         bias_diff_complete_cu = bias_cu - bias_cu_complete)

if (all(is.na(metrics$rel_efficiency_cu))) {
  warning("rel_efficiency_cu is NA everywhere - is the complete_data arm collected?")
}

saveRDS(metrics, file.path(study$res_path, "bin_miss_metrics.RDS"))

# BLP and independence tests run on the true CATE instead of an estimated one
# (true nuisances too - see run_true_cate_tests() in R/cate_models.R), to see
# how the tests themselves perform independent of any estimator's error.
# multiple_imputation rows come back NA/NA - see true_cate_test_row()'s doc
# and README.md's "still has no HTE tests, for any model" note.
true_cate_tests <- compute_run_metrics(study, all_results_df, true_cate_test_row)
saveRDS(true_cate_tests, file.path(study$res_path, "bin_miss_true_cate_tests.RDS"))

print("metrics calculated!")
