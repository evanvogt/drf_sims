##########
# title: metrics for binary outcome confidence intervals
##########
# The pipeline lives in R/metrics.R. Grouping columns come from the study
# config, so this script carries only what is specific to binary outcome confidence intervals.

library(here)
source(here("confidence_intervals/binary/bin_ci_config.R"))
source(here("R", "metrics.R"))

all_results_df <- readRDS(file.path(study$res_path, "ci_bin_all.RDS"))

metrics <- compute_metrics(
  study, all_results_df, models = CI_MODELS,
  per_model = function(model_res, true_tau, model, sim_res, keys) {
    out <- bind_cols(tibble(model = model),
                     interval_metrics(model_res$hb_lb, model_res$hb_ub, true_tau))

    # the causal forest also reports its own variance, so score that interval
    # too - as a separate row, labelled causal_forest_inbuilt
    if (model == "causal_forest" && !is.null(model_res$variance)) {
      ni <- normal_interval(model_res$tau, model_res$variance)
      out <- bind_rows(out, bind_cols(tibble(model = "causal_forest_inbuilt"),
                                      interval_metrics(ni$lb, ni$ub, true_tau)))
    }

    # grid-based (as opposed to per-unit) simultaneous coverage, when the
    # replicate carries one - see R/dgm_scenarios.R::build_query_grid and
    # bin_ci_analysis.R. Own row, labelled "<model>_grid", the same pattern
    # causal_forest_inbuilt uses to distinguish an alternative interval for
    # the same underlying arm.
    if (!is.null(model_res$grid_lb)) {
      out <- bind_rows(out, bind_cols(
        tibble(model = paste0(model, "_grid")),
        interval_metrics(model_res$grid_lb, model_res$grid_ub, sim_res$grid_truth$tau)
      ))
    }

    out
  }
)

saveRDS(metrics, file.path(study$res_path, "ci_bin_metrics.RDS"))
print("metrics calculated!")
