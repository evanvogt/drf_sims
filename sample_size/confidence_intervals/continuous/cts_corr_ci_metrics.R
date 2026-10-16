##########
# title: metrics for correlated continuous outcome confidence intervals
##########
# As confidence_intervals/continuous/cts_ci_metrics.R. rho is one of the
# study's path_cols, so every metrics row carries it as a key column, and the
# bias-eliminated grid reference is averaged within each rho too.

library(here)
source(here("sample_size/correlated/confidence_intervals/continuous/cts_corr_ci_config.R"))
source(here("R", "metrics.R"))

all_results_df <- readRDS(file.path(study$res_path, "cts_corr_ci_all.RDS"))

# Bias-eliminated (BE) reference: across-run mean of tau_grid at each fixed
# grid point, per (rho, scenario, n, CI_sf, model) cell - see
# R/metrics.R::grid_be_reference().
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

    # grid-based (as opposed to per-unit) simultaneous coverage - own row,
    # labelled "<model>_grid", as the parent study
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

saveRDS(metrics, file.path(study$res_path, "cts_corr_ci_metrics.RDS"))
print("metrics calculated!")
