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
# three result trees without an edit. With 8 score types (see me_score_types()
# below) x 2 pipelines x however many arms:
#
#   main        cv, whole                                    -> 1 + 8x2x2 = 33
#   strategies  whole, cv_indep, cv_shared, holdout           -> 1 + 8x4x2 = 65
#   split       split                                         -> 1 + 8x1x2 = 17
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

# ---- the propensity axis ----------------------------------------------------
# Every score below comes in two versions: one using the nuisance pipeline's
# ESTIMATED propensity (df$pi, df$phi) and one with the propensity FIXED at
# 0.5 (the literal 0.5 passed to calc_infl_score, df$phi05). The fixed version
# is not an approximation here - R/dgm_scenarios.R assigns treatment with
# `W <- rbinom(n, 1, 0.5)`, a fair coin independent of X, in every scenario of
# every set. So 0.5 is the TRUE propensity, the fixed-pi scores are oracle-pi
# scores, and `dr` vs `dr05` isolates exactly one thing: the cost of estimating
# a propensity that never needed estimating.
#
# Nothing had to be rerun to add them. phi05 has been in
# me_nuisance.R::calculate_pseudos() since this study's first commit, so every
# completed replicate already carries it, and calc_infl_score() already took pi
# as an argument. See README.md's "Scores" section.

#' Influence-corrected PEHE proxy - unchanged from metric_utils.R
#'
#' @param pi estimated propensity, or the scalar 0.5 for the oracle-pi twin.
#'   Every operation here is vectorised arithmetic, so a length-1 pi recycles
#'   to give exactly the same answer as rep(0.5, n) - me_testing.R checks this
#'   rather than leaving it to be assumed.
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

#' DR calibration score
#'
#' Splits the evaluation rows into k_groups quantile groups of the CANDIDATE's
#' own tau_hat, and asks whether the group average treatment effect the
#' candidate implies matches the one the DR scores imply:
#'
#'   GATE_k^hat = mean of tau_hat over G_k     (what the candidate claims)
#'   GATE_k^DR  = mean of phi over G_k         (what the pseudo-outcomes say)
#'   M          = sum_k |G_k| * |GATE_k^hat - GATE_k^DR|
#'
#' A loss, lower is better, like every other column - so me_results.qmd's
#' ranking/regret machinery needs no special case for it.
#'
#' TWO DELIBERATE DEPARTURES from the obvious formula, both easy to "fix" back:
#'
#'  1. ABSOLUTE, not signed. The signed sum sum_k |G_k| (GATE_k^hat - GATE_k^DR)
#'     telescopes to sum_i (tau_hat_i - phi_i) = n (mean(tau_hat) - mean(phi)):
#'     every group boundary cancels and what is left is an ATE-bias measure that
#'     cannot see miscalibration at all. Over- and under-estimation in different
#'     groups must not be allowed to net out, which is what the |.| is for.
#'  2. WEIGHTS ARE GROUP COUNTS, NOT PROPORTIONS - sum(w) is n, not 1. That
#'     makes the score scale with n. Within one (scenario, n, run) the factor is
#'     the same for all 9 candidates, so it moves no ranking, no Spearman
#'     correlation against true_pehe, no pick and no regret; what it does mean is
#'     that the raw score is not comparable ACROSS sample sizes. See README.md.
#'
#' WHAT THIS SCORE CANNOT SEE. Unlike calc_infl_score()/calc_dr_risk() it is not
#' an estimate of PEHE, and it measures calibration rather than discrimination.
#' A candidate predicting the constant ATE for every unit has arbitrary quantile
#' groups, so each group's GATE^DR is about the overall mean phi and every
#' discrepancy is near zero - a near-perfect score at whatever PEHE the true
#' heterogeneity implies. Expect the shrunk-to-a-constant candidates (net2, and
#' everything in scenario 1) to look good here. That is a property of the
#' metric, not a bug to patch, and it is part of what this study is measuring.
#'
#' @param tau_hat the candidate's CATE predictions on the evaluation rows
#' @param phi the DR score on the same rows - df$phi, or df$phi05 for the
#'   oracle-pi twin
#' @param k_groups number of quantile groups (see CAL_QUANTILES in me_config.R)
calc_cal_score <- function(tau_hat, phi, k_groups) {
  n <- length(tau_hat)

  # Grouped by RANK, not by cut(tau_hat, quantile(tau_hat, ...)). Cut points
  # collapse to fewer than k_groups groups - or emit NA rows - as soon as
  # tau_hat is constant or heavily tied, and this study produces exactly that:
  # scenario 1 has no heterogeneity, and a fully-shrunk net* candidate predicts
  # one number for every row. Ranking with ties broken by position always gives
  # k_groups non-empty, near-equal groups, so the score stays defined and the
  # column never silently becomes NA for the flattest candidates - the ones
  # whose calibration is most worth reading.
  g <- ceiling(k_groups * rank(tau_hat, ties.method = "first") / n)

  gate_hat <- tapply(tau_hat, g, mean)
  gate_dr <- tapply(phi, g, mean)
  w <- as.numeric(table(g)) # counts, not proportions - see note 2 above

  sum(w * abs(gate_hat - gate_dr))
}

#' The score types one (pipeline, arm) contributes, in column order
#'
#' The `05` in a name means the propensity was fixed at 0.5; `q<K>` is the
#' calibration score's group count. Both are folded into this single token
#' rather than added as further underscore-separated fields, because
#' me_results.qmd recovers the design axes from the column name as
#' <score_type>_<arm>_<pipeline> - a fourth field would break that parse. No
#' score type contains "05" unless its propensity is fixed (calq10 has a "10",
#' never an "05"), so the pi regime is recoverable with grepl("05", .).
#'
#' Exists so me_testing.R can assert the column count against the same source
#' me_per_model() builds from, instead of hardcoding a number that drifts.
me_score_types <- function(k_groups = CAL_QUANTILES) {
  c(
    "infl", "infl05",
    "dr", "dr05",
    paste0("calq", k_groups),
    paste0("cal05q", k_groups)
  )
}

#' Every score one (pipeline, arm) contributes for one candidate
#'
#' Names must match me_score_types() exactly and in order; me_testing.R checks
#' that they do, which is what keeps the two definitions from drifting apart.
arm_scores <- function(tau_hat, df, Y, W, k_groups = CAL_QUANTILES) {
  # A loud failure beats scoring an arm whose phi05 was never computed and
  # silently emitting garbage. Every arm goes through
  # me_nuisance.R::calculate_pseudos(), which has always produced this column.
  stopifnot(
    "nuisance data.frame has no phi05 - see me_nuisance.R::calculate_pseudos()" =
      "phi05" %in% names(df)
  )

  c(
    list(
      infl   = calc_infl_score(tau_hat, df$tau_T, df$pi, Y, W),
      infl05 = calc_infl_score(tau_hat, df$tau_T, 0.5, Y, W),
      dr     = calc_dr_risk(tau_hat, df$phi),
      dr05   = calc_dr_risk(tau_hat, df$phi05)
    ),
    setNames(
      lapply(k_groups, function(k) calc_cal_score(tau_hat, df$phi, k)),
      paste0("calq", k_groups)
    ),
    setNames(
      lapply(k_groups, function(k) calc_cal_score(tau_hat, df$phi05, k)),
      paste0("cal05q", k_groups)
    )
  )
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
      arm <- arm_scores(tau_hat, df, Y, W)
      names(arm) <- paste0(names(arm), "_", fold_type, "_", pipeline)
      scores <- c(scores, arm)
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
