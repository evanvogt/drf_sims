##########
# title: metrics for missing covariates, continuous outcome
##########
# The pipeline lives in R/metrics.R. Grouping columns come from the study
# config, so this script carries only what is specific to missing covariates, continuous outcome.

library(here)
source(here("missing/continuous/cts_miss_config.R"))
source(here("R", "metrics.R"))

all_results_df <- readRDS(file.path(study$res_path, "cts_miss_all.RDS"))

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

saveRDS(metrics, file.path(study$res_path, "cts_miss_metrics.RDS"))
print("metrics calculated!")
