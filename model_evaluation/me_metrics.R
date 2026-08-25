##########
# title: metrics for the model evaluation study
##########
# Reuses R/metrics.R::compute_metrics() rather than a bespoke aggregator.
# This works because me_analysis.R saves the 9 candidate models as TOP-LEVEL
# keys of each run's results (rf1..SL3), exactly like continuous saves
# causal_forest/dr_random_forest/etc. - so compute_metrics()'s
# intersect(names(sim_res), models) finds them directly, and its per_model()
# callback receives the full per-run result, so it can reach
# sim_res$nuisances/sim_res$data while scoring one candidate at a time.
#
# The proxy-score formulas (calc_infl_score, calc_dr_risk) are unchanged from
# the old metric_utils.R - only the orchestration moved, from one vectorized
# pass over all 9 models inside calc_metrics(), to compute_metrics()'s own
# per-(run, model) loop calling me_per_model() once per model. Same numbers;
# me_testing.R's regression check confirms this.
#
# extract_surrogate()/calc_scores_for_config()/fix_automl() from the old
# metric_utils.R are gone: fix_automl() guarded against a "p1"-named 5th
# column that predict_automl_models() (me_nuisance.R) never actually
# produces - confirmed by reading it directly, it always names that column
# "pi" - and extract_surrogate()'s cv/infold/whole-wide surrogate matrix is
# unnecessary once scoring pulls df$tau_T/df$pi/df$phi directly from one
# fold-type's data.frame at a time.
#
# The score columns are DERIVED from whatever arms the nuisance list carries,
# never enumerated here - which is what lets the same me_per_model() serve all
# three result trees without an edit:
#
#   main        cv, whole                                    -> 1 + 2x2x2 =  9
#   strategies  whole, cv_indep, cv_shared, holdout           -> 1 + 2x4x2 = 17
#   split       split                                         -> 1 + 2x1x2 =  5
#
# The split tree works unchanged for a second reason too: me_split.R stores
# `data` and `truth` ALREADY RESTRICTED to its 20% evaluation rows, so
# tau_hat, Y, W, the nuisance data.frame and true_tau are all the same length
# there, and every formula below is length-agnostic. Scoring an 80:20 arm on
# its holdout therefore needed no per_model variant, only a results object
# assembled with the restriction already applied.
#
# Note: compute_metrics() always does true_tau <- sim_res$truth$tau - there
# is no equivalent of the old calc_metrics(..., truth_avail = FALSE) branch.
# me_analysis.R only ever runs against dgm_scenarios.R-generated data (which
# always carries truth), so this is a real but low-risk scope reduction.

library(here)
source(here("model_evaluation/me_config.R"))
source(here("R", "metrics.R"))

#' Influence-corrected PEHE proxy - unchanged from metric_utils.R
calc_infl_score <- function(tau_hat, tau, pi, Y, W) {
  A <- W - pi
  C <- pi * (1 - pi)
  B <- 2 * W * A / C
  diff_tau <- tau_hat - tau

  pehe_infl <- (1 - B) * tau^2 + B * Y * diff_tau - A * diff_tau^2 + tau_hat^2

  mean(pehe_infl)
}

#' DR risk proxy (MSE against the AIPW pseudo-outcome) - unchanged from metric_utils.R
calc_dr_risk <- function(tau_hat, tau_dr) {
  mean((tau_hat - tau_dr)^2)
}

#' Every proxy score, plus true PEHE, for one candidate model of one run
#'
#' Matches the R/metrics.R::compute_metrics() per_model() contract:
#' function(model_res, true_tau, model, sim_res, keys).
me_per_model <- function(model_res, true_tau, model, sim_res, keys) {
  tau_hat <- model_res$tau
  Y <- sim_res$data$Y
  W <- sim_res$data$W

  scores <- list(true_pehe = mean((tau_hat - true_tau)^2))

  for (pipeline in names(sim_res$nuisances)) { # xgb, automl
    # whatever arms this tree carries - see the table above. Deliberately not
    # checked against NUISANCE_ARMS: the same function has to score the main
    # tree's cv/whole and the split tree's single arm too.
    for (fold_type in names(sim_res$nuisances[[pipeline]])) {
      df <- sim_res$nuisances[[pipeline]][[fold_type]]
      scores[[paste0("infl_", fold_type, "_", pipeline)]] <- calc_infl_score(
        tau_hat, df$tau_T, df$pi, Y, W
      )
      scores[[paste0("dr_", fold_type, "_", pipeline)]] <- calc_dr_risk(
        tau_hat, df$phi
      )
    }
  }

  tibble::as_tibble(scores)
}

# Which tree to score:
#
#   Rscript me_metrics.R                # the main study
#   Rscript me_metrics.R strategies     # the nuisance-arm pass
#   Rscript me_metrics.R split          # the 80:20 arm
#
# Unrecognised arguments fall back to "main" rather than erroring, because
# me_testing.R sources this file for its function definitions and would
# otherwise hand its own `full` argument to me_study().
which_study <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(which_study) || !which_study %in% c("main", "strategies", "split")) {
  which_study <- "main"
}
st <- me_study(which_study)

# The compute-and-save step only runs once there's something to collect. This
# guard lets me_testing.R source this file purely for calc_infl_score()/
# calc_dr_risk()/me_per_model() without needing the collected file to exist.
all_path <- file.path(st$res_path, paste0(st$prefix, "_all.RDS"))
if (file.exists(all_path)) {
  all_results_df <- readRDS(all_path)

  metrics <- compute_metrics(
    st, all_results_df, models = CANDIDATE_MODELS, per_model = me_per_model
  )

  saveRDS(metrics, file.path(st$res_path, paste0(st$prefix, "_metrics.RDS")))
  print(paste0("metrics calculated! (", st$name, ")"))
} else {
  message(
    basename(all_path), " not found at ", all_path,
    " - sourced for its function definitions only. Run me_collect.R first, ",
    "then `Rscript me_metrics.R` directly, to actually compute metrics."
  )
}
