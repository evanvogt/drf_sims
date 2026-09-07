##########
# title: metrics for continuous outcome confidence intervals
##########
# The pipeline lives in R/metrics.R. Grouping columns come from the study
# config, so this script carries only what is specific to continuous outcome confidence intervals.

library(here)
source(here("confidence_intervals/continuous/cts_ci_config.R"))
source(here("R", "metrics.R"))

all_results_df <- readRDS(file.path(study$res_path, "ci_cts_all.RDS"))

# Bias-eliminated (BE) reference: across-run mean of tau_grid at each fixed
# grid point, per (scenario, n, CI_sf, model) cell - substitutes for the
# unknown true tau in BE-coverage. Grid-only: see R/metrics.R::grid_be_reference()
# for why the per-unit hb_lb/hb_ub arm has no analogous reference. Computed
# once and captured by the per_model() closure below; compute_metrics() itself
# is unchanged.
be_ref <- grid_be_reference(study, all_results_df, models = CI_MODELS)

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
    # cts_ci_analysis.R. Own row, labelled "<model>_grid", the same pattern
    # causal_forest_inbuilt uses to distinguish an alternative interval for
    # the same underlying arm.
    if (!is.null(model_res$grid_lb)) {
      out <- bind_rows(out, bind_cols(
        tibble(model = paste0(model, "_grid")),
        interval_metrics(model_res$grid_lb, model_res$grid_ub, sim_res$grid_truth$tau),
        be_interval_metrics(model_res$grid_lb, model_res$grid_ub,
                            be_reference_for(be_ref, keys, model, study$path_cols))
      ))
    }

    out
  }
)

saveRDS(metrics, file.path(study$res_path, "ci_cts_metrics.RDS"))
print("metrics calculated!")
