##########
# title: shared machinery for the interim-analysis validation studies
##########
# Everything the continuous/ and binary/ arms do identically: splitting one
# trial at the interim, fitting the three estimators and both importance
# measures on a chunk, the TE-VIM / TreeSHAP / interaction-test helpers, and the
# four chunk-1-vs-chunk-2 comparisons. Each arm's <prefix>_val_models.R sources
# this and fixes the outcome family; its <prefix>_val_analysis.R is the array
# entry point.
#
# This lived in continuous/cts_val_models.R and continuous/cts_val_analysis.R
# until the binary arm was added (2026-10-08) and moved here unchanged, apart
# from two arguments that default to the continuous behaviour:
#   family  the DR SuperLearner's outcome-model family (fit_val_methods)
#   robust  HC3 standard errors in the interaction tests (interaction_pval,
#           interaction_pval_adj, chunk_validations) - see interaction_pval
#
# The estimation logic here (both importance measures and the interaction
# tests) is specific to these studies - nothing else in the repo uses it - so it
# stays in validation/ rather than moving into R/.
#
# The nuisance estimation and the three estimators (causal_forest,
# dr_random_forest, dr_superlearner) live in R/cate_models.R. The question is
# whether subgroups/variance/variable-importance found on one chunk of a trial
# replicate on the next, not which estimator is best, so the oracle and
# semi-oracle arms are not run.
#
# Adopting R/cate_models.R's run_causal_forest also switches the causal-forest
# arm to the shared Y.hat.cf nuisance (a regression of Y on X alone) rather than
# the old validation-local nuisance (Y on (W,X) with the observed W plugged
# back in). That is a deliberate change, not an accident - see
# continuous/README.md's Status section.
#
# The two forest estimators (and their TE-VIM refits below) are whole-sample
# OOB, not fold-crossfit - see R/cate_models.R's crossfitting-strategy note.
# n_folds is used by the DR SuperLearner alone: its nuisances, its stage 2 and
# its TE-VIM refits all share one set of folds.

library(xgboost)
library(SHAPforxgboost)
library(rpart)

source(here::here("R", "cate_models.R"))

###################
# The trial and its chunks
###################

#' Split one generated trial at the interim analysis
#'
#' The first n * interim_prop participants are chunk 1, the rest chunk 2. Rows
#' are iid, so the first n1 are a valid "enrolled by the interim" cohort. The
#' chunks used to be drawn as two separate datasets, and
#' generate_scenario_data() calibrates bW to 80% power at the n it is given - so
#' each chunk was its own trial with its own ATE, not two halves of one.
#' Splitting one draw also nests the chunks across interim_prop: run r's chunk 1
#' at 0.25 is a subset of its chunk 1 at 0.30, so curves over interim_prop are
#' paired within run.
#'
#' round() because n * interim_prop is not always an exact integer in floating
#' point (1000 * 0.35 need not be 350), and an index sequence would truncate it.
#'
#' @param gen generate_scenario_data()'s output for the whole trial
#' @return list(data1, data2, truth1, truth2), chunk 2's row names reset
split_trial <- function(gen, interim_prop) {
  n <- nrow(gen$dataset)
  n1 <- round(n * interim_prop)
  chunk1 <- seq_len(n1)
  chunk2 <- (n1 + 1):n

  data2 <- gen$dataset[chunk2, ]
  rownames(data2) <- NULL
  truth2 <- gen$truth[chunk2, , drop = FALSE]
  rownames(truth2) <- NULL

  list(data1 = gen$dataset[chunk1, ], data2 = data2,
       truth1 = gen$truth[chunk1, , drop = FALSE], truth2 = truth2)
}

#' DR SuperLearner folds for a chunk of m rows
#'
#' The forests are whole-sample OOB and use none. 5 folds below 500 rows, 10
#' from 500; chunks run from 250 to 750 rows, so a fold's held-out set is never
#' under 50.
chunk_folds <- function(m) if (m < 500) 5L else 10L

###################
# Fitting one chunk
###################

#' Fit the three estimators and both importance measures on one trial chunk
#'
#' @param family the DR SuperLearner's outcome-model family - gaussian() for a
#'   continuous outcome, binomial() for a binary one (SuperLearner then uses
#'   method.NNloglik). The forests take a binary Y as it is: their outcome
#'   regressions estimate the risk directly, as in R/cate_models.R.
#' @param num.threads grf thread count, forwarded to every regression_forest()/
#'   causal_forest() call below - including the TE-VIM refits, which run inside
#'   `future_map()` and so would otherwise each grab every visible core. NULL
#'   keeps grf's own default (all cores).
#' @param verbose_timing if TRUE, time each block below with R/utils.R's
#'   `timed()` and attach the elapsed seconds as `results$timings`. Same idea as
#'   R/cate_models.R::cate_methods(), and the same caveat: `timings` is then an
#'   element of `results` that is not a model, so callers deriving a model list
#'   from `names(results)` must exclude it (chunk_validations() does).
#' @param sl_lib the DR SuperLearner's libraries, list(W = , Y = , tau = ).
#'   Defaults to R/sl_library.R's sl_libraries() at this chunk's size, which is
#'   what the analysis scripts pass too.
fit_val_methods <- function(data, n_folds = 10, num.threads = NULL,
                            verbose_timing = FALSE,
                            sl_lib = sl_libraries(nrow(data)),
                            family = gaussian()) {

  X <- as.matrix(data[, -c(1:2)])
  Y <- data$Y
  W <- data$W

  # contiguous folds, as R/cate_models.R::cate_methods() builds them
  fold_indices <- sort(seq(nrow(X)) %% n_folds) + 1
  fold_list <- unique(fold_indices)

  timings <- list()
  time_step <- function(label, expr) {
    if (verbose_timing) {
      t <- timed(expr)
      timings[[label]] <<- t$time
      t$value
    } else {
      expr
    }
  }

  cat("Computing nuisance functions...\n")
  nuisances <- time_step("nuisance_rf", nuisance_rf(X, Y, W, num.threads = num.threads))

  results <- list()

  cat("Running Causal Forest...\n")
  results$causal_forest <- time_step(
    "causal_forest", run_causal_forest(X, Y, W, nuisances, num.threads = num.threads))
  results$causal_forest$te_vims <- time_step("cf_te_vims", get_te_vims_causal_forest(
    X, Y, W, nuisances$po, results$causal_forest$tau, num.threads = num.threads
  ))
  results$causal_forest$shap_vims <- time_step(
    "cf_shap_vims", get_shap_vims(X, results$causal_forest$tau))

  cat("Running DR Random Forest...\n")
  results$dr_random_forest <- time_step(
    "dr_random_forest",
    run_dr_random_forest(X, Y, W, nuisances, num.threads = num.threads))
  results$dr_random_forest$te_vims <- time_step("dr_te_vims", get_te_vims(
    X, nuisances$po, results$dr_random_forest$tau, num.threads = num.threads
  ))
  results$dr_random_forest$shap_vims <- time_step("dr_shap_vims", get_shap_vims(
    X, results$dr_random_forest$tau
  ))

  # SuperLearner wants a data.frame; the forests above keep the matrix. Its
  # pseudo-outcome comes from its own SuperLearner nuisances, not nuisance_rf's,
  # so its TE-VIMs are scored against that.
  cat("Running DR SuperLearner...\n")
  X_df <- as.data.frame(X)
  nuisances_sl <- time_step("nuisance_sl",
                            nuisance_sl(X_df, Y, W, fold_indices, sl_lib,
                                        family = family))
  dropped_nuis <- attr(nuisances_sl, "sl_dropped")
  attr(nuisances_sl, "sl_dropped") <- NULL

  results$dr_superlearner <- time_step("dr_superlearner", run_dr_superlearner(
    X_df, Y, W, nuisances_sl, fold_indices, fold_list, sl_lib))
  results$dr_superlearner$te_vims <- time_step("sl_te_vims", get_te_vims_superlearner(
    X_df, nuisances_sl$po, results$dr_superlearner$tau, fold_indices, sl_lib,
    results$dr_superlearner$sl_dropped
  ))
  results$dr_superlearner$shap_vims <- time_step("sl_shap_vims", get_shap_vims(
    X, results$dr_superlearner$tau
  ))
  # kept inside the model's own element: every top-level element of `results`
  # other than data/truth/timings is read as a model by chunk_validations()
  results$dr_superlearner$sl_dropped <- rbind(dropped_nuis,
                                              results$dr_superlearner$sl_dropped)

  if (verbose_timing) results$timings <- timings

  results
}

###################
# TE-VIMs (treatment-effect variable importance measures)
###################
# For each covariate, refit the whole-sample stage-2 model with that
# covariate dropped and take its OOB predictions, then compare their
# prediction error against the full model's OOB predictions - a
# pathwise-derivative-style importance measure. Whole-sample OOB throughout,
# matching R/cate_models.R's stage2_whole_rf / run_causal_forest - this used
# to refit per fold instead, before those moved off fold-crossfitting.
# `get_te_vims` covers estimators whose second stage is a plain regression
# forest on a pseudo-outcome (DR Random Forest); `get_te_vims_causal_forest`
# refits a causal_forest instead, since that estimator takes (X, Y, W)
# directly rather than (X, pseudo-outcome). `get_te_vims_superlearner` is the
# fold-crossfit exception: the DR SuperLearner's stage 2 is crossfit, so its
# dropped-covariate refits are too, on the same folds.

#' Score dropped-covariate CATE predictions against the full model's
#'
#' The one TE-VIM calculation every refit function below shares: the increase
#' in squared error against the pseudo-outcome from dropping each covariate,
#' with its influence-function standard error.
#'
#' @param po pseudo-outcome vector
#' @param tau full-model out-of-sample CATE estimates
#' @param sub_taus n x p matrix, column j the predictions without covariate j,
#'   columns named by covariate
#' @return data.frame, rows tevim and std_err, one column per covariate
te_vim_scores <- function(po, tau, sub_taus) {
  n_obs <- length(po)
  r_tau <- (po - tau)^2

  te_vims <- apply(sub_taus, 2, function(sub_tau) {
    r_subtau <- (po - sub_tau)^2
    tevim <- sum(r_subtau - r_tau) / n_obs
    infl <- r_subtau - r_tau - tevim
    std_err <- sqrt(sum(infl^2)) / n_obs
    list(tevim = tevim, std_err = std_err)
  }) %>% simplify2array()

  as.data.frame(te_vims)
}

#' TE-VIMs for a DR-learner style estimator (pseudo-outcome regression)
#'
#' @param X covariate matrix
#' @param po pseudo-outcome vector
#' @param tau full-model OOB CATE estimates
#' @param num.threads grf thread count for each dropped-covariate refit. These
#'   run one per worker, so leaving it NULL means every worker spawns a forest on
#'   all visible cores at once, oversubscribing whatever PBS allocated.
get_te_vims <- function(X, po, tau, num.threads = NULL) {
  covariates <- colnames(X)

  sub_taus_list <- future_map(seq_along(covariates), function(i) {
    new_X <- as.matrix(X[, -i])
    predict(regression_forest(new_X, po, num.threads = num.threads))$predictions
  }, .options = furrr_options(seed = TRUE))

  sub_taus <- do.call(cbind, sub_taus_list)
  colnames(sub_taus) <- covariates

  te_vim_scores(po, tau, sub_taus)
}

#' TE-VIMs for the causal forest, refitting a whole-sample causal_forest per
#' dropped covariate
#'
#' Each dropped-covariate forest cross-fits its own nuisances internally, same
#' as the full model (R/cate_models.R::run_causal_forest) - this used to
#' refit fold-wise against externally-supplied double-crossfit nuisances.
#'
#' @param po the po field of nuisance_rf()'s output, used for scoring only -
#'   the forest itself no longer takes externally-supplied nuisances
#' @param tau full-model OOB CATE estimates
#' @param num.threads grf thread count, as in get_te_vims() - and it matters more
#'   here, since a causal_forest cross-fits its own nuisances and so is the more
#'   expensive of the two refits
get_te_vims_causal_forest <- function(X, Y, W, po, tau, num.threads = NULL) {
  covariates <- colnames(X)

  sub_taus_list <- future_map(seq_along(covariates), function(i) {
    new_X <- as.matrix(X[, -i])
    predict(causal_forest(new_X, Y, W, num.threads = num.threads))$predictions
  }, .options = furrr_options(seed = TRUE))

  sub_taus <- do.call(cbind, sub_taus_list)
  colnames(sub_taus) <- covariates

  te_vim_scores(po, tau, sub_taus)
}

#' TE-VIMs for the DR SuperLearner, refitting its crossfit stage 2 per dropped
#' covariate
#'
#' Same pseudo-outcome, same folds, same stage-2 library as the full model
#' (R/cate_models.R::stage_2_sl), with covariate j removed - so tau and every
#' sub_tau are out-of-fold predictions of the same learner. Each fold reuses the
#' library the full model's pretest left on that fold rather than pretesting
#' again: the pretest only removes learners that error or predict non-finite
#' values, which dropping one covariate does not change, and skipping it saves
#' 9 two-fold fits per refit. The p x fold refits are one flat future_map, so
#' extra workers spread across all of them. The stage 2 is a regression of the
#' pseudo-outcome, so it is gaussian() whatever the outcome.
#'
#' @param X covariate data.frame, as SuperLearner takes it
#' @param po nuisance_sl()'s pseudo-outcome
#' @param tau the full model's crossfit CATE estimates
#' @param fold_indices the folds the full model used
#' @param sl_lib list(W = , Y = , tau = ); the tau library is used
#' @param sl_dropped run_dr_superlearner()'s dropped-learner table, whose
#'   model == "tau" rows say what the full model's pretest removed per fold
get_te_vims_superlearner <- function(X, po, tau, fold_indices, sl_lib,
                                     sl_dropped = NULL) {
  covariates <- colnames(X)
  tau_lib <- as_sl_libs(sl_lib)$tau

  fold_lib <- function(fold) {
    gone <- if (is.null(sl_dropped)) character() else
      sl_dropped$learner[sl_dropped$model == "tau" & sl_dropped$fold == fold]
    lib <- setdiff(tau_lib, gone)
    if (length(lib) == 0) "SL.mean" else lib
  }

  jobs <- expand.grid(j = seq_along(covariates), fold = unique(fold_indices))

  preds <- future_map(seq_len(nrow(jobs)), function(k) {
    j <- jobs$j[k]
    fold <- jobs$fold[k]
    in_train <- fold_indices != fold
    new_X <- X[, -j, drop = FALSE]
    fit <- sl_fit_predict(po[in_train], new_X[in_train, , drop = FALSE],
                          list(tau = new_X[!in_train, , drop = FALSE]),
                          fold_lib(fold), family = gaussian())
    fit$pred$tau
  }, .options = furrr_options(seed = TRUE))

  sub_taus <- matrix(NA_real_, length(po), length(covariates),
                     dimnames = list(NULL, covariates))
  for (k in seq_len(nrow(jobs))) {
    sub_taus[fold_indices == jobs$fold[k], jobs$j[k]] <- preds[[k]]
  }

  te_vim_scores(po, tau, sub_taus)
}


###################
# Surrogate TreeSHAP
###################
# The second variable-importance measure. Where the TE-VIMs above ask "how much
# worse does the second stage predict without this covariate", TreeSHAP asks
# "how much of this unit's estimated CATE is attributable to this covariate",
# and averages |that| over units. The two need not agree, which is the point of
# carrying both: chunk_validations() compares each one's chunk-1 vs chunk-2
# ranking separately.
#
# This is "Strategy 3" / indirect SHAP from Svensson et al.'s SHAP_CATE - fit a
# tuned xgboost surrogate to the *estimated* CATEs, then take exact TreeSHAP of
# the surrogate. Because it reads only (X, tau) it applies unchanged to every
# estimator, unlike the TE-VIMs which need an estimator-specific refit. The
# CATE is continuous for either outcome (a risk difference for a binary one), so
# the surrogate is a squared-error regression in both arms.
#
# Adapted from SHAP_CATE-main/FUNCTIONS/functions.R::cvboost3(), with four
# deliberate departures, all for runtime - this runs six times per array job
# across 1100 jobs:
#   - 3 learning rates x 3 depths (9 combos), not 8 x 3 (24)
#   - 5-fold CV and at most 2000 rounds, not 10-fold and 10000
#   - the tuning grid is parallelised over the workers the analysis driver has
#     already planned, with nthread = 1 inside each fit
#   - the xgb.cv prediction callback is dropped. Only the final xgb.train model
#     is needed for SHAP and the round count comes off the evaluation log, so
#     saving the fold models was pure memory - and it also avoids the callback
#     having been renamed cb.cv.predict -> xgb.cb.cv.predict at xgboost 2.0,
#     which would otherwise tie this to a version nothing in this repo pins.

CVBOOST_ETAS <- c(0.01, 0.05, 0.1)
CVBOOST_DEPTHS <- 2:4

#' Cross-validated xgboost surrogate for a continuous target
#'
#' @param x feature matrix
#' @param y target vector (here, the estimated CATEs)
#' @return the fitted xgb.Booster at the best (eta, max_depth, nrounds)
cvboost_cate <- function(x, y, k_folds = 5, ETAs = CVBOOST_ETAS,
                         TR.DEPTH = CVBOOST_DEPTHS, ntrees_max = 2000,
                         early_stopping_rounds = 50) {
  tune_grid <- expand.grid(eta = ETAs, max_depth = TR.DEPTH)

  base_param <- function(eta, max_depth) {
    list(objective = "reg:squarederror", eval_metric = "rmse",
         subsample = 0.9, colsample_bytree = 0.9,
         eta = eta, max_depth = max_depth, nthread = 1)
  }

  fits <- future_map(seq_len(nrow(tune_grid)), function(i) {
    param <- base_param(tune_grid$eta[i], tune_grid$max_depth[i])
    cvfit <- xgboost::xgb.cv(
      data = xgboost::xgb.DMatrix(data = x, label = y),
      params = param,
      nfold = k_folds,
      nrounds = ntrees_max,
      early_stopping_rounds = early_stopping_rounds,
      maximize = FALSE,
      verbose = FALSE
    )
    # which.min rather than $best_iteration: same answer, no dependence on a
    # field whose presence has moved around between xgboost versions
    losses <- as.data.frame(cvfit$evaluation_log)[["test_rmse_mean"]]
    list(param = param, loss = min(losses), nrounds = which.min(losses))
  }, .options = furrr_options(seed = TRUE))

  best <- fits[[which.min(vapply(fits, function(f) f$loss, numeric(1)))]]

  xgboost::xgb.train(
    data = xgboost::xgb.DMatrix(data = x, label = y),
    params = best$param,
    nrounds = best$nrounds
  )
}

#' Surrogate-TreeSHAP variable importance for a fitted CATE
#'
#' Shaped to match get_te_vims()'s output so chunk_validations() can read either
#' with the same `[1, ]` indexing: one column per covariate, in `colnames(X)`
#' order. Only the mean-|SHAP| summary is returned, not the n x p per-unit SHAP
#' matrix - nothing downstream uses the latter and it would add ~200KB per run
#' to <prefix>_val_all.RDS.
#'
#' @param X covariate matrix
#' @param tau estimated CATEs to explain
get_shap_vims <- function(X, tau) {
  surrogate <- cvboost_cate(X, tau)

  shap <- SHAPforxgboost::shap.values(xgb_model = surrogate, X_train = X)

  # shap.values() returns mean_shap_score sorted descending by importance, not
  # in X's column order. Reindexing is not cosmetic: chunk_validations() ranks
  # this against te_vims by position, so leaving it sorted would silently
  # scramble every rank.
  out <- as.data.frame(as.list(shap$mean_shap_score[colnames(X)]))
  colnames(out) <- colnames(X)
  rownames(out) <- "mean_abs_shap"
  out
}


###################
# Interaction tests
###################

#' p-value of one coefficient of a fitted lm, classical or HC3
#'
#' @param robust FALSE: summary.lm's t-test, as the continuous arm has always
#'   used. TRUE: the same t-test on sandwich::vcovHC(type = "HC3") standard
#'   errors. A binary outcome makes `Y ~ W * v` a linear probability model,
#'   whose errors are heteroskedastic by construction (Var = p(1 - p)), so the
#'   classical standard errors are wrong there.
#' @return the p-value, or NA_real_ if lm dropped the term (aliased) or its
#'   robust standard error is not finite (a leverage-1 row under HC3)
coef_pval <- function(fit, term, robust = FALSE) {
  if (!robust) {
    co <- summary(fit)$coefficients
    if (!term %in% rownames(co)) return(NA_real_)
    return(unname(co[term, 4]))
  }
  est <- coef(fit)
  if (!term %in% names(est) || is.na(est[[term]])) return(NA_real_)
  V <- sandwich::vcovHC(fit, type = "HC3")
  se <- sqrt(V[term, term])
  if (!is.finite(se) || se <= 0) return(NA_real_)
  unname(2 * stats::pt(-abs(est[[term]] / se), df = fit$df.residual))
}

#' Interaction p-value for `Y ~ W * v` in the later chunk
#'
#' Indexes the coefficient by name, not position. The positional form is what
#' made bottom_pval report the *intercept* p-value (`pvals_bottom[1]`) rather
#' than the W-by-subgroup interaction, and it is fragile for a second reason:
#' when `v` is constant, lm drops the interaction and the coefficient matrix has
#' fewer than four rows, so `[4]` reads whatever happens to be there.
#'
#' For a binary outcome the model is linear in the risk, which is the scale the
#' binary DGMs put the treatment effect on (R/dgm_scenarios.R, BINARY OUTCOMES
#' ARE ON THE RISK-DIFFERENCE SCALE), so W:v is a risk-difference interaction.
#'
#' @param v subgroup indicator or covariate to interact with treatment
#' @param robust HC3 standard errors - see coef_pval()
#' @return the W:v p-value, or NA_real_ if `v` carries no contrast (e.g. rpart
#'   predicted no bottom10 leaf into chunk 2)
interaction_pval <- function(Y, W, v, robust = FALSE) {
  d <- data.frame(Y = Y, W = W, v = v)
  if (length(unique(na.omit(d$v))) < 2) return(NA_real_)
  coef_pval(lm(Y ~ W * v, data = d), "W:v", robust)
}

#' Interaction p-value for one covariate, adjusted for every other covariate
#'
#' The W:x_top coefficient of `Y ~ W * (all covariates)` - does x_top modify the
#' effect once the other covariates' main effects and interactions are held
#' fixed? interaction_pval() above is the marginal test, and under correlated
#' covariates the two answer different questions: with rho = 0.5, E[tau | X5]
#' moves with X5 because X5 is correlated with X4, so a marginal W x X5 test
#' finds an interaction even though X5 modifies nothing. A wrong top covariate
#' then "replicates" marginally; this adjusted test is the one that should not.
#'
#' @param X data frame of every covariate, x_top among its columns
#' @param x_top name of the covariate whose interaction is tested
#' @param robust HC3 standard errors - see coef_pval()
#' @return the W:x_top p-value, or NA_real_ if lm dropped that term (aliased)
interaction_pval_adj <- function(Y, W, X, x_top, robust = FALSE) {
  d <- data.frame(Y = Y, W = W, X)
  coef_pval(lm(Y ~ W * ., data = d), paste0("W:", x_top), robust)
}


###################
# The four chunk comparisons
###################

#' Compare what chunk 1 found against chunk 2
#'
#' @param results1,results2 fit_val_methods() on each chunk. Every element other
#'   than data/truth/timings is read as a model.
#' @param data1,data2 the chunks themselves
#' @param robust HC3 standard errors in every interaction test - see coef_pval()
#' @return list(subgroups, variances, var_imps, top_var_tests), each a named
#'   list by model
chunk_validations <- function(results1, results2, data1, data2, robust = FALSE) {

  X1 <- data1[, -c(1, 2)]
  X2 <- data2[, -c(1, 2)]

  # "timings" is only present when fit_val_methods() was called with
  # verbose_timing = TRUE, which the analysis scripts do not do, but it is
  # excluded anyway so every non-data element of `results1` is a model.
  models <- setdiff(names(results1), c("data", "truth", "timings"))

  ##########
  # subgroups based on top and bottom responders
  ##########
  # Fit a tree on chunk 1's top/bottom-10%-CATE groups, predict them into chunk
  # 2, and test the W x subgroup interaction there. interaction_pval indexes the
  # W:v coefficient by name. This used to be two positional lookups, and the
  # bottom one read `pvals_bottom[1]` - the intercept - so bottom_pval was never
  # a subgroup test at all. Results from before that fix are not comparable.
  subgroups <- list()
  for (model in models) {
    tau1 <- results1[[model]]$tau

    group <- cut(tau1,
                 breaks = quantile(tau1, probs = c(0, 0.1, 0.9, 1)),
                 labels = c("bottom10", "middle", "top10"),
                 include.lowest = TRUE)
    df_train <- data.frame(group = group, X1)

    tree_group <- rpart(group ~ ., data = df_train, method = "class")
    group_pred <- predict(tree_group, newdata = X2, type = "class")

    subgroups[[model]] <- c(
      top    = interaction_pval(data2$Y, data2$W, as.numeric(group_pred == "top10"),
                                robust),
      bottom = interaction_pval(data2$Y, data2$W, as.numeric(group_pred == "bottom10"),
                                robust)
    )
  }

  ##########
  # Compare variance between early and later chunks
  ##########
  variances <- list()
  for (model in models) {
    vt1 <- var(results1[[model]]$tau)
    vt2 <- var(results2[[model]]$tau)
    variances[[model]] <- c(vt1 = unname(vt1), vt2 = unname(vt2))
  }

  ##########
  # Compare variable importance between early and late chunks
  ##########
  # Two measures: the TE-VIMs and surrogate TreeSHAP. Both are
  # larger-is-more-important, so rank() means the same thing for each - rank 1
  # is the least important covariate.
  measure_fields <- c(tevim = "te_vims", shap = "shap_vims")

  var_imps <- list()
  for (model in models) {
    fit1 <- results1[[model]]
    fit2 <- results2[[model]]

    per_measure <- lapply(names(measure_fields), function(measure) {
      field <- measure_fields[[measure]]
      imp1 <- unlist(fit1[[field]][1, ])
      imp2 <- unlist(fit2[[field]][1, ])

      data.frame(variables = colnames(fit1[[field]]),
                 measure = measure,
                 vi1 = rank(imp1),
                 vi2 = rank(imp2),
                 stringsAsFactors = FALSE) %>%
        mutate(diff = vi2 - vi1)
    })

    var_imps[[model]] <- do.call(rbind, per_measure)
  }

  ##########
  # Carry the top-ranked covariate into the remaining participants
  ##########
  # The point of ranking covariates is whether the winner means anything, so
  # take each measure's chunk-1 most important covariate and interaction-test it
  # in chunk 2: the continuous W x X_top interaction (works for continuous and
  # already-binary covariates alike, no arbitrary cut point), the same adjusted
  # for every other covariate (p_cts_adj - with correlated covariates a
  # non-modifier correlated with X4 shows a real *marginal* interaction), and a
  # median split, parallel to the top10/bottom10 tests above. x_top2 is chunk
  # 2's own winner, kept so the report can ask how often the two chunks even
  # agree on which covariate matters most.
  top_var_tests <- list()
  for (model in models) {
    vi <- var_imps[[model]]

    rows <- lapply(split(vi, vi$measure), function(v) {
      x_top <- v$variables[which.max(v$vi1)]
      xt <- data2[[x_top]]

      data.frame(measure = v$measure[1],
                 x_top = x_top,
                 x_top2 = v$variables[which.max(v$vi2)],
                 p_cts = interaction_pval(data2$Y, data2$W, xt, robust),
                 p_cts_adj = interaction_pval_adj(data2$Y, data2$W, X2, x_top, robust),
                 p_split = interaction_pval(data2$Y, data2$W,
                                            as.numeric(xt > median(xt)), robust),
                 stringsAsFactors = FALSE)
    })

    top_var_tests[[model]] <- do.call(rbind, rows)
  }

  # TODO: all three estimators carry BLP_whole/independence_cate/independence_po
  # (see R/cate_models.R, R/metrics.R::hte_test_metrics()) in a shape a
  # chunk-vs-chunk HTE-test comparison could use directly. Not implemented yet -
  # see continuous/README.md's Status section.

  list(subgroups = subgroups, variances = variances,
       var_imps = var_imps, top_var_tests = top_var_tests)
}
