##########
# title: metrics for correlated binary outcome
##########
# As binary/bin_metrics.R. rho is one of the study's path_cols, so every
# metrics row carries it as a key column alongside scenario and n.

library(here)
source(here("sample_size/correlated/binary/bin_corr_config.R"))
source(here("R", "metrics.R"))
source(here("R", "cate_models.R"))

all_results_df <- readRDS(file.path(study$res_path, "bin_corr_all.RDS"))

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

saveRDS(metrics, file.path(study$res_path, "bin_corr_metrics.RDS"))

# BLP and independence tests run on the true CATE instead of an estimated one
# (true nuisances too - see run_true_cate_tests() in R/cate_models.R)
true_cate_tests <- compute_run_metrics(study, all_results_df, true_cate_test_row)
saveRDS(true_cate_tests, file.path(study$res_path, "bin_corr_true_cate_tests.RDS"))

print("metrics calculated!")
