##########
# title: shared metric definitions
##########
# Definitions only, no side effects: each study's *_metrics.R sources this and
# supplies its own grouping columns.
#
# Taken from crossfitting/cf_metrics.R, which was the best-written of the six
# near-identical copies. The others differed only in which grouping columns they
# carried and in guarding against absent HTE tests.

require(dplyr)
require(tibble)
require(tidyr)   # unnest_longer
require(purrr)   # map, map_int

# ---- bias sign --------------------------------------------------------------
# TEMPORARY (bug G). Historically `bias` was mean(true - est) while `ate_bias`
# in the same tibble was mean(est) - mean(true), so the two disagreed in sign.
# The agreed fix is `est - true` everywhere. Until the metrics are regenerated,
# the default reproduces the old behaviour so the refactor can be proved inert;
# Step 8 flips this to "est-true" and then deletes the option entirely.
BIAS_SIGN <- "est-true"

#' Point metrics for one set of CATE estimates against the known truth
#'
#' @param est estimated CATEs
#' @param true true CATEs
#' @param scenario scenario index. Scenario 1 has no heterogeneity, so the true
#'   CATE is constant, the correlation metrics are undefined, and both are
#'   reported as 0 rather than NA - the convention the study has always used.
#' @param bias_sign "est-true" (the convention) or "true-est" (legacy). Affects
#'   `bias` only; `ate_bias` has always been est - true.
cate_metrics <- function(est, true, scenario, bias_sign = BIAS_SIGN) {
  bias_sign <- match.arg(bias_sign, c("true-est", "est-true"))
  bias <- if (bias_sign == "est-true") {
    mean(est - true, na.rm = TRUE)
  } else {
    mean(true - est, na.rm = TRUE)
  }

  tibble(
    bias = bias,
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

#' Heterogeneity-test p-values recorded alongside a fitted model
#'
#' Guarded, because the CI studies skip the tests and the missing-data studies
#' skip them whenever the covariate matrix still has NAs.
#'
#' @param model_res one model's entry in a per-run results object
hte_test_metrics <- function(model_res) {
  tibble(
    BLP_p = if (!is.null(model_res$BLP_whole)) {
      model_res$BLP_whole[4, 2]
    } else NA_real_,
    indep_cate = if (!is.null(model_res$independence_cate)) {
      as.numeric(model_res$independence_cate$p_value)
    } else NA_real_,
    indep_po = if (!is.null(model_res$independence_po)) {
      as.numeric(model_res$independence_po$p_value)
    } else NA_real_
  )
}

#' Every metric for one model of one run
#'
#' @param model_res one model's entry in a per-run results object
#' @param true true CATEs for that run
#' @param scenario scenario index
run_model_metrics <- function(model_res, true, scenario, bias_sign = BIAS_SIGN) {
  bind_cols(
    cate_metrics(model_res$tau, true, scenario, bias_sign),
    hte_test_metrics(model_res)
  )
}

# the CATE estimators the point metrics apply to, in display order
CATE_MODELS <- c("causal_forest", "dr_random_forest", "dr_oracle",
                 "dr_semi_oracle", "dr_superlearner")

# the CI studies drop the SuperLearner arm
CI_MODELS <- c("causal_forest", "dr_random_forest", "dr_oracle", "dr_semi_oracle")

#' Coverage and width of one interval estimate
#'
#' Marginal coverage is the proportion of units their own interval covers;
#' simultaneous coverage is whether the band covers every unit at once, which is
#' what the half-sample bootstrap is constructed to control.
interval_metrics <- function(lb, ub, true) {
  tibble(
    marginal_coverage = mean(as.numeric(true >= lb & true <= ub)),
    simultaneous_coverage = as.numeric(all(true >= lb & true <= ub)),
    mean_ci_length = mean(ub - lb)
  )
}

#' Normal-approximation interval from a variance estimate
#'
#' Used for the causal forest's own variance estimates, alongside the bootstrap.
normal_interval <- function(tau, variance, alpha = 0.05) {
  se <- sqrt(variance)
  list(lb = tau + qnorm(alpha / 2) * se,
       ub = tau + qnorm(1 - alpha / 2) * se)
}

#' Turn a collected results tibble into one metrics row per (combination, run, model)
#'
#' Replaces seven near-identical unnest / pmap / map_dfr pipelines. The grouping
#' columns are taken from the study config rather than retyped, so a study that
#' gains a factor does not need its metrics script edited.
#'
#' @param study the study config (supplies path_cols)
#' @param all_results_df output of get_results()
#' @param models which model names to look for inside each run
#' @param per_model function(model_res, true_tau, model, sim_res, keys) returning
#'   a tibble of that model's metric columns. `keys` is the one-row data frame of
#'   grouping values, so `keys$scenario` is available to cate_metrics(). The
#'   callback may return SEVERAL rows and may include its own `model` column,
#'   which is how the CI studies add the "causal_forest_inbuilt" row alongside
#'   the bootstrap one.
compute_metrics <- function(study, all_results_df, models = CATE_MODELS,
                            per_model) {

  group_cols <- study$path_cols

  df <- all_results_df %>%
    unnest_longer(results) %>%
    mutate(run = map_int(results, ~ .x$run),
           sim_res = map(results, ~ .x$result)) %>%
    select(-results)

  keys <- df[, c(group_cols, "run"), drop = FALSE]

  rows <- lapply(seq_len(nrow(df)), function(i) {
    sim_res <- df$sim_res[[i]]
    true_tau <- sim_res$truth$tau
    models_run <- intersect(names(sim_res), models)

    key_row <- keys[i, , drop = FALSE]

    bind_rows(lapply(models_run, function(m) {
      out <- per_model(sim_res[[m]], true_tau, m, sim_res, key_row)
      if (!"model" %in% names(out)) out <- bind_cols(tibble(model = m), out)
      bind_cols(key_row[rep(1, nrow(out)), , drop = FALSE], out)
    }))
  })

  bind_rows(rows)
}
