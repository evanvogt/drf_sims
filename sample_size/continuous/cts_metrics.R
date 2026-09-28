##########
# title: metrics for continuous outcome
##########
# The pipeline lives in R/metrics.R. Grouping columns come from the study
# config, so this script carries only what is specific to continuous outcome.

library(here)
source(here("sample_size/continuous/cts_config.R"))
source(here("R", "metrics.R"))
source(here("R", "cate_models.R"))

all_results_df <- readRDS(file.path(study$res_path, "cts_all.RDS"))

# The T-learners are derived here from each run's saved DR stage 1, not fitted
# at estimation time - see add_t_learners() in R/metrics.R.
metrics <- compute_metrics(
  study, all_results_df, models = c(CATE_MODELS, T_LEARNER_MODELS),
  augment = add_t_learners,
  per_model = function(model_res, true_tau, model, sim_res, keys) {
    bind_cols(
      cate_metrics(model_res$tau, true_tau, keys$scenario),
      hte_test_metrics(model_res, sim_res, model)
    )
  }
)

saveRDS(metrics, file.path(study$res_path, "cts_metrics.RDS"))

# BLP and independence tests run on the true CATE instead of an estimated one
# (true nuisances too - see run_true_cate_tests() in R/cate_models.R), to see
# how the tests themselves perform independent of any estimator's error.
true_cate_tests <- compute_run_metrics(study, all_results_df, true_cate_test_row)
saveRDS(true_cate_tests, file.path(study$res_path, "cts_true_cate_tests.RDS"))

print("metrics calculated!")
