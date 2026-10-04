##########
# title: metrics for the correlated optimal sample.fraction calibration (continuous)
##########
# Scores each run's final band - built at the run's own calibrated CI_sf -
# against the true CATE, exactly as the CI studies score theirs
# (R/metrics.R::interval_metrics, be_interval_metrics). Nominal coverage is
# 0.90 here (alpha = 0.1 in the analysis script), not the CI studies' 0.95.
#
# A run saves its single arm (dr_random_forest) at the top level rather than
# under a model name, so as_model_result() below files it under
# sim_res$dr_random_forest first. That lets the shared compute_metrics() and
# grid_be_reference() run unchanged.

library(here)
source(here("sample_size/correlated/confidence_intervals/optimal_sf/cts_corr_ci_sf_config.R"))
source(here("R", "metrics.R"))

SF_MODEL <- "dr_random_forest"

#' File a calibration run's top-level estimate under sim_res$dr_random_forest
as_model_result <- function(sim_res) {
  sim_res[[SF_MODEL]] <- sim_res[c("tau", "hb_lb", "hb_ub",
                                   "tau_grid", "grid_lb", "grid_ub")]
  sim_res
}

#' The calibration's own view of the pick: mean coverage of tau.hat (the
#' plug-in target) and mean width at the chosen CI_sf
pick_columns <- function(sim_res) {
  cc <- sim_res$coverage_curve
  k <- which.min(abs(cc$sf - sim_res$optimal_sf))
  tibble(optimal_sf = sim_res$optimal_sf,
         plugin_coverage = cc$mean_coverage[k],
         plugin_ci_width = cc$mean_ci_width[k])
}

all_results_df <- readRDS(file.path(study$res_path, "cts_corr_ci_sf_all.RDS"))
all_results_df$results <- lapply(all_results_df$results, function(runs) {
  lapply(runs, function(r) {
    r$result <- as_model_result(r$result)
    r
  })
})

# Bias-eliminated (BE) reference: across-run mean of tau_grid at each fixed
# grid point, per (rho, scenario, n) cell - see R/metrics.R::grid_be_reference()
be_ref <- grid_be_reference(study, all_results_df, models = SF_MODEL)

metrics <- compute_metrics(
  study, all_results_df, models = SF_MODEL,
  per_model = function(model_res, true_tau, model, sim_res, keys) {
    out <- bind_cols(tibble(model = model),
                     interval_metrics(model_res$hb_lb, model_res$hb_ub, true_tau))

    if (!is.null(model_res$grid_lb)) {
      out <- bind_rows(out, bind_cols(
        tibble(model = paste0(model, "_grid")),
        interval_metrics(model_res$grid_lb, model_res$grid_ub, sim_res$grid_truth$tau),
        be_interval_metrics(model_res$grid_lb, model_res$grid_ub,
                            be_reference_for(be_ref, keys, model, study$path_cols))
      ))
    }

    bind_cols(out, pick_columns(sim_res)[rep(1, nrow(out)), ])
  }
)

saveRDS(metrics, file.path(study$res_path, "cts_corr_ci_sf_metrics.RDS"))
print("metrics calculated!")
