##########
# title: metrics for binary outcome
##########
# The pipeline lives in R/metrics.R. Grouping columns come from the study
# config, so this script carries only what is specific to binary outcome.

library(here)
source(here("binary/bin_config.R"))
source(here("R", "metrics.R"))

all_results_df <- readRDS(file.path(study$res_path, "bin_all.RDS"))

metrics <- compute_metrics(
  study, all_results_df, models = CATE_MODELS,
  per_model = function(model_res, true_tau, model, sim_res, keys) {
    bind_cols(
      cate_metrics(model_res$tau, true_tau, keys$scenario),
      hte_test_metrics(model_res)
    )
  }
)

saveRDS(metrics, file.path(study$res_path, "bin_metrics.RDS"))
print("metrics calculated!")
