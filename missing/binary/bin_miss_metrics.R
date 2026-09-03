##########
# title: metrics for missing covariates, binary outcome
##########
# The pipeline lives in R/metrics.R. Grouping columns come from the study
# config, so this script carries only what is specific to missing covariates, binary outcome.

library(here)
source(here("missing/binary/bin_miss_config.R"))
source(here("R", "metrics.R"))
source(here("R", "cate_models.R"))

all_results_df <- readRDS(file.path(study$res_path, "bin_miss_all.RDS"))

metrics <- compute_metrics(
  study, all_results_df, models = CATE_MODELS,
  per_model = function(model_res, true_tau, model, sim_res, keys) {
    bind_cols(
      cate_metrics(model_res$tau, true_tau, keys$scenario),
      hte_test_metrics(model_res)
    )
  }
)

# Relative efficiency and relative bias against the complete-data reference arm.
# This used to be NA for every row (bug C): the reference was selected with
# method == "complete_data", which the collect grid did not contain. It does now.
# rel_ate_bias and rel_bias_cate (relative to the TRUE parameter, not this
# complete_data arm) already arrive from cate_metrics() in R/metrics.R - no
# join needed for those.
ref_by <- setdiff(c(study$path_cols, "run", "model"), "method")

complete_ref <- metrics %>%
  filter(method == "complete_data") %>%
  select(all_of(ref_by), mse_complete = mse, bias_complete = bias)

metrics <- metrics %>%
  left_join(complete_ref, by = ref_by) %>%
  mutate(rel_efficiency = mse / mse_complete,
         rel_bias_complete = bias / bias_complete)

if (all(is.na(metrics$rel_efficiency))) {
  warning("rel_efficiency is NA everywhere - is the complete_data arm collected?")
}

if (all(is.na(metrics$rel_bias_complete))) {
  warning("rel_bias_complete is NA everywhere - is the complete_data arm collected?")
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
