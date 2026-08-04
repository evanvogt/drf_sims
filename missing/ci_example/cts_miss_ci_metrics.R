##########
# title: metrics for MI confidence-interval example
##########
# The pipeline lives in R/metrics.R. Grouping columns come from the study
# config, so this script carries only what is specific to MI confidence-interval example.

library(here)
source(here("missing/ci_example/cts_miss_ci_config.R"))
source(here("R", "metrics.R"))

all_results_df <- readRDS(file.path(study$res_path, "cts_miss_ci_all.RDS"))

# Three MI pooling strategies are compared, so each model contributes three
# rows, labelled by strategy. See combine_mi_ci() in R/bootstrap_ci.R.
STRATEGIES <- c(pooled = "pooled", mib = "mib", hybrid = "hybrid")

metrics <- compute_metrics(
  study, all_results_df, models = CI_MODELS,
  per_model = function(model_res, true_tau, model, sim_res, keys) {
    bind_rows(lapply(names(STRATEGIES), function(s) {
      lb <- model_res[[paste0("lb_", s)]]
      ub <- model_res[[paste0("ub_", s)]]
      if (is.null(lb) || is.null(ub)) return(NULL)
      bind_cols(tibble(model = model, strategy = s),
                interval_metrics(lb, ub, true_tau),
                cate_metrics(model_res$tau, true_tau, keys$scenario))
    }))
  }
)

saveRDS(metrics, file.path(study$res_path, "cts_miss_ci_metrics.RDS"))
print("metrics calculated!")
