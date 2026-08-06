##########
# title: metric definitions for the crossfitting comparison
##########
# Definitions only - no side effects. cf_collect.R sources this and applies it
# while streaming the per-run files, so the metrics tibble is built in one pass
# and the large nested intermediate the other studies write is never materialised.
#
# Note on sign: bias is (estimate - truth) here, the usual convention. cts_metrics.R:37
# uses mean(true - est) for bias while its ate_bias at line 43 uses the opposite
# sign, so these numbers are not sign-comparable with the main continuous study.

require(dplyr)
require(tibble)

#' Point metrics for one set of CATE estimates against the known truth
#'
#' @param est estimated CATEs
#' @param true true CATEs
#' @param scenario scenario index; scenario 1 has no heterogeneity so the
#'   correlation metrics are undefined and forced to 0, as in cts_metrics.R:41
cate_metrics <- function(est, true, scenario) {
  tibble(
    bias = mean(est - true, na.rm = TRUE),
    ate_bias = mean(est, na.rm = TRUE) - mean(true, na.rm = TRUE),
    mse = mean((est - true)^2, na.rm = TRUE),
    rmse = sqrt(mean((est - true)^2, na.rm = TRUE)),
    mae = mean(abs(est - true), na.rm = TRUE),
    corr = if (scenario != 1) cor(true, est, use = "pairwise.complete.obs") else 0,
    spearman = if (scenario != 1) {
      cor(true, est, method = "spearman", use = "pairwise.complete.obs")
    } else 0,
    sign_acc = mean(sign(est) == sign(true), na.rm = TRUE),
    n_na = sum(is.na(est))
  )
}

#' Metrics for every arm of one replicate, on both the training and test samples
#'
#' @param sim_res one per-run object as written by cf_analysis.R
#' @param scenario scenario index
#' @return tibble, one row per (arm, set)
run_metrics <- function(sim_res, scenario) {
  bind_rows(lapply(names(sim_res$arms), function(nm) {
    a <- sim_res$arms[[nm]]

    train <- cate_metrics(a$tau, sim_res$truth_tau, scenario)
    test <- cate_metrics(a$tau_test, sim_res$truth_test_tau, scenario)

    # `%||%` is only base R from 4.4.0 and the cluster runs 4.3.2
    single <- if (is.null(a$mse_test_single)) NA_real_ else a$mse_test_single

    bind_rows(mutate(train, set = "train"), mutate(test, set = "test")) %>%
      mutate(
        scenario = scenario,
        run = sim_res$run,
        arm = nm,
        family = a$family,
        variant = a$variant,
        time_nuisance = a$time_nuisance,
        time_stage2 = a$time_stage2,
        time_total = a$time_nuisance + a$time_stage2,
        # like-for-like test score: each fold model scored separately then
        # averaged, rather than the V-model ensemble that mse scores. equals mse
        # on the test row for whole-sample arms. see arm() in cf_models.R
        mse_test_single = single,
        .before = 1
      )
  }))
}

# display order and labels, used by cf_results.R
variant_levels <- c("dcf", "scf_scf", "scf_scf_new", "scf_oob",
                    "scf_oob_t", "oob_oob", "oob_oob_s", "oob_oob_manual",
                    "cf_dcf", "cf_scf", "cf_full_oob", "cf_default")

variant_labels <- c(
  dcf = "Double CF",
  scf_scf = "CF+CF",
  scf_scf_new = "CF+CF (new)",
  scf_oob = "CF+OOB",
  scf_oob_t = "CF+OOB (T)",
  oob_oob = "OOB+OOB",
  oob_oob_s = "OOB+OOB (S)",
  oob_oob_manual = "OOB+OOB (S, manual)",
  cf_dcf = "Double CF nuis.",
  cf_scf = "CF nuis.",
  cf_full_oob = "Full+OOB",
  cf_default = "grf OOB"
)

family_labels <- c(dr_rf = "DR-learner (RF)",
                   dr_sl = "DR-learner (SuperLearner)",
                   causal_forest = "Causal forest")
