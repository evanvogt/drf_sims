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

#' Point metrics for one set of CATE estimates against the known truth
#'
#' @param est estimated CATEs
#' @param true true CATEs
#' @param scenario scenario index. Scenario 1 has no heterogeneity, so the true
#'   CATE is constant, the correlation metrics are undefined, and both are
#'   reported as 0 rather than NA - the convention the study has always used.
#'
#' `rel_ate_bias` and `rel_bias_cate` are "relative bias" in the conventional,
#' vs-true-parameter sense (as opposed to `rel_efficiency`/`rel_bias_complete`
#' in the missing-data studies, which are ratios against the complete-data
#' arm). `rel_ate_bias` is the parameter-level version (`ate_bias` divided by
#' the true ATE); `rel_bias_cate` is the per-unit version, `(est - true) /
#' true` averaged over units. True CATE is heterogeneous and crosses zero in
#' several scenarios, so both guard a zero denominator with NA rather than
#' letting it produce Inf/NaN.
cate_metrics <- function(est, true, scenario) {
  bias <- mean(est - true, na.rm = TRUE)
  ate_bias <- mean(est, na.rm = TRUE) - mean(true, na.rm = TRUE)
  true_ate <- mean(true, na.rm = TRUE)
  cate_ratio <- (est - true) / true
  cate_ratio[true == 0] <- NA_real_

  tibble(
    bias = bias,
    ate_bias = ate_bias,
    rel_ate_bias = if (true_ate == 0) NA_real_ else ate_bias / true_ate,
    rel_bias_cate = mean(cate_ratio, na.rm = TRUE),
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
run_model_metrics <- function(model_res, true, scenario) {
  bind_cols(
    cate_metrics(model_res$tau, true, scenario),
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

#' Bias-eliminated coverage of an interval against a substituted reference
#'
#' Same construction as interval_metrics()'s coverage columns, but scored
#' against a stand-in reference (e.g. the across-run mean of the point
#' estimate at a fixed query point) instead of the unknown true parameter -
#' "BE-coverage": https://joonho112.github.io/simsum-mini-course/06-metrics-inference.html#sec-becoverage
#' Substituting the mean estimate for the truth isolates whether the
#' interval's width is correctly calibrated, independent of whether the point
#' estimator itself is biased.
#'
#' No internal na.rm, matching interval_metrics() - one NA in `be_reference`
#' or the bounds makes the whole row NA. mean_ci_length is deliberately not
#' recomputed here: lb/ub are whatever the caller already scored with
#' interval_metrics(), so it would be a byte-for-byte duplicate column -
#' callers bind_cols() this alongside that call's result on the same row.
#'
#' @param lb,ub the interval bounds, as passed to interval_metrics()
#' @param be_reference the bias-eliminated reference vector (same length as
#'   lb/ub), or NULL when none is available (e.g. this model/cell never
#'   produced one) - returns NA columns rather than erroring, the same
#'   graceful-absence pattern hte_test_metrics() uses for missing test objects.
be_interval_metrics <- function(lb, ub, be_reference) {
  if (is.null(be_reference)) {
    return(tibble(be_marginal_coverage = NA_real_,
                  be_simultaneous_coverage = NA_real_))
  }
  tibble(
    be_marginal_coverage = mean(as.numeric(be_reference >= lb & be_reference <= ub)),
    be_simultaneous_coverage = as.numeric(all(be_reference >= lb & be_reference <= ub))
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

#' One row per (run, sim_res) with this study's grouping keys attached
#'
#' Factored out of compute_metrics() so a driver that doesn't need a per-model
#' loop (e.g. compute_run_metrics()) can share the same unnesting.
#'
#' @param study the study config (supplies path_cols)
#' @param all_results_df output of get_results()
#' @return list(df = one row per run with a sim_res list-column, keys = the
#'   matching grouping-column + run data frame, row-aligned to df)
unnest_results <- function(study, all_results_df) {
  df <- all_results_df %>%
    unnest_longer(results) %>%
    mutate(run = map_int(results, ~ .x$run),
           sim_res = map(results, ~ .x$result)) %>%
    select(-results)

  list(df = df, keys = df[, c(study$path_cols, "run"), drop = FALSE])
}

#' Bias-eliminated grid reference: across-run mean of tau_grid at each fixed
#' grid point, per (path_cols..., model) cell
#'
#' The theta_bar substitute BE-coverage uses in place of the unknown true
#' theta - see be_interval_metrics(). Meaningful only where the grid rows
#' apply (confidence_intervals/{continuous,binary}): the query grid
#' (R/dgm_scenarios.R::build_query_grid()) is fixed within a (scenario, n)
#' pair, and calibrate_bW() is a deterministic function of (params, n) that
#' consumes no RNG, so every run of one path_cols cell targets the exact same
#' grid_truth$tau - unlike the per-unit hb_lb/hb_ub arm, where a fresh sample
#' of units is drawn every run and there is no fixed per-unit truth to
#' average toward. Grouped at the full path_cols granularity (not pooled
#' across CI_sf, even though CI_sf never enters bW/Z_query/grid_truth)
#' because CI_sf does enter the half-sample bootstraps run between arms
#' within one replicate, so a later arm's tau_grid cannot be assumed
#' identical across CI_sf for the same run.
#'
#' `tau_grid` is present for a model in a run iff `grid_lb` is (both come
#' from the same `!is.null(Z_query)` gate in R/cate_models.R::cate_methods()
#' and R/bootstrap_ci.R's rf_oob_half_boot()/cf_oob_half_boot()), so
#' filtering on tau_grid here exactly matches the `!is.null(model_res$grid_lb)`
#' gate the *_ci_metrics.R scripts already use to decide whether to emit a
#' "<model>_grid" row.
#'
#' @param study the study config (supplies path_cols)
#' @param all_results_df output of get_results()
#' @param models which model names to average tau_grid for
#' @return named list keyed by "<path_cols pasted with \\r>\\r<model>", each
#'   value a numeric vector the same length as that cell's Z_query. A
#'   (cell, model) that never produced a tau_grid is simply absent from the
#'   list - see be_reference_for().
grid_be_reference <- function(study, all_results_df, models = CI_MODELS) {
  u <- unnest_results(study, all_results_df)
  cell_key <- do.call(paste, c(u$keys[, study$path_cols, drop = FALSE], sep = "\r"))

  acc <- list()
  for (i in seq_len(nrow(u$df))) {
    sim_res <- u$df$sim_res[[i]]
    for (m in intersect(names(sim_res), models)) {
      tg <- sim_res[[m]]$tau_grid
      if (is.null(tg)) next
      k <- paste(cell_key[i], m, sep = "\r")
      acc[[k]] <- c(acc[[k]], list(tg))
    }
  }

  lapply(acc, function(vecs) rowMeans(do.call(cbind, vecs), na.rm = TRUE))
}

#' Look up one run's bias-eliminated reference vector
#'
#' Keys the same way grid_be_reference() built its lookup, so a per_model()
#' closure that already has `keys` (the one-row grouping data frame
#' compute_metrics() hands it) and `model` can fetch its cell's reference
#' without recomputing anything.
#'
#' @param be_ref output of grid_be_reference()
#' @param keys one-row keys data frame, as compute_metrics() passes to per_model()
#' @param model model name
#' @param path_cols study$path_cols
#' @return numeric vector, or NULL if this (cell, model) never contributed to
#'   be_ref (see grid_be_reference())
be_reference_for <- function(be_ref, keys, model, path_cols) {
  k <- paste(do.call(paste, c(keys[1, path_cols, drop = FALSE], sep = "\r")),
             model, sep = "\r")
  be_ref[[k]]
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

  u <- unnest_results(study, all_results_df)
  df <- u$df
  keys <- u$keys

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

#' One row per run from a callback that only needs the run's own sim_res
#'
#' The per-run counterpart of compute_metrics(): for metrics that aren't tied
#' to any one fitted model (e.g. a test on the true CATE), skipping the
#' per-model loop entirely.
#'
#' @param study the study config (supplies path_cols)
#' @param all_results_df output of get_results()
#' @param test_fn function(sim_res) returning a one-row tibble of that run's
#'   metric columns
compute_run_metrics <- function(study, all_results_df, test_fn) {
  u <- unnest_results(study, all_results_df)

  rows <- lapply(seq_len(nrow(u$df)), function(i) {
    bind_cols(u$keys[i, , drop = FALSE], test_fn(u$df$sim_res[[i]]))
  })

  bind_rows(rows)
}
