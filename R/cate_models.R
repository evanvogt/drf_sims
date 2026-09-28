##########
# title: shared CATE estimators (DR-learner family + causal forest)
##########
# One implementation of what used to be seven near-identical copies:
#   sample_size/continuous/cts_models.R              sample_size/binary/bin_models.R
#   missing/continuous/cts_miss_models.R missing/binary/bin_miss_models.R
#   missing/ci_example/cts_miss_ci_models.R
#   sample_size/confidence_intervals/continuous/cts_ci_models.R
#   sample_size/confidence_intervals/binary/bin_ci_models.R
#
# Those copies differed on four axes. Three are the arguments below:
#
#   family   gaussian() vs binomial(), controlling the SuperLearner outcome model
#            (family + method.NNloglik).
#   ipw      sample.weights (grf) / obsWeights (SuperLearner) for the missing-data
#            IPW arm. NULL reproduces the unweighted path exactly.
#   ci       list(boot=, sf=, alpha=) turns on the half-sample bootstrap.
#   profile  which historical variant's *orchestration* to reproduce - see below.
#
# The fourth, where the inverse link lives for the oracle arm, is gone. Every
# oracle formula in R/dgm_scenarios.R returns the outcome MEAN - the linear
# predictor for a continuous outcome, the risk for a binary one - so the oracle
# arm applies no link. Until the binary DGM moved to the risk-difference scale
# the binary formulas were linear predictors, and an `oracle_link` argument said
# whether to apply plogis; missing/binary/ once passed "identity" to formulas
# that no longer carried their plogis (bug M).
#
# `profile` exists because the variants also disagree about which post-estimation
# tests get run, and those disagreements look like drift rather than design. They
# are reproduced exactly rather than harmonised, because harmonising would change
# results for studies that are not otherwise re-running. See PROFILES.
#
# Crossfitting strategy, per crossfitting/'s comparison of alternatives against
# the double-crossfitting this file used to do throughout:
#   dr_random_forest, dr_oracle, dr_semi_oracle  whole-sample OOB, "oob_oob":
#     a T-learner outcome model (t_learner_rf - one forest per arm) with no
#     sample splitting, each unit's own-arm prediction out-of-bag and its
#     other-arm prediction from a forest that never saw it, and an OOB stage-2
#     regression forest (stage2_whole_rf). Until the move to per-arm outcome
#     models this was "oob_oob_s", an S-learner forest on cbind(W, X) read at
#     counterfactual rows through grf's X.orig (oob_predict_counterfactual).
#   causal_forest   grf's own internal cross-fitting, "cf_default": a plain
#     causal_forest(X, Y, W) with no externally-supplied nuisances.
#   dr_superlearner   single leave-one-fold-out crossfit, "scf_scf": nuisance_sl
#     and stage_2_sl share the same fold_indices, rather than double-crossfit
#     nuisances feeding a separately-split stage 2. Its outcome model is also a
#     T-learner: one SuperLearner per arm (libraries in R/sl_library.R).
# See crossfitting/cf_models.R for the full arm comparison these three replace.

require(coin)
require(grf)
require(SuperLearner)
require(GenericML)
require(future)
require(furrr)
require(dplyr)

source(here::here("R", "utils.R"))        # collate_predictions, setup_rng_stream
source(here::here("R", "bootstrap_ci.R")) # cf_half_boot, rf_half_boot
source(here::here("R", "sl_library.R"))   # as_sl_libs, sl_fit_predict, pretest_superlearner

# ---- flags pending a decision ----------------------------------------------


# ---- orchestration profiles -------------------------------------------------
#
#                        base      ci        missing
#  causal forest variance  no       yes       yes
#  causal forest tests     yes      no        yes
#  dr_random_forest tests  yes      no        yes     <- was "no", see NOTE
#  oracle / semi tests     yes      no        yes
#  SuperLearner arm        yes      no        yes (skipped if X still has NAs)
#  half-sample bootstrap   no       yes       no
#
# NOTE: the missing profile used to set dr_rf_tests = FALSE, so that arm alone
# was built inline as list(tau = stage2_whole_rf(...)$tau) and carried no BLP or
# independence test while every other arm did - copy-paste drift from the file
# missing/ was forked from, not a decision. The decision was then taken that
# every model should carry the tests where possible, and this is it.
#
# Results made before the change were back-filled in place by a one-off patch
# script (removed after commit e7b1d59; the results it patched are archived,
# and patched files lack dr_random_forest$variance). The re-run needs no patch.
#
# The CI profiles keep tests off; that one IS deliberate (see
# sample_size/confidence_intervals/README.md).
#
# aggregate_nuisances (row-mean summaries of a double-crossfit matrix) used to
# be a fourth axis here; it was dropped when nuisance_rf/nuisance_sl moved to
# whole-sample OOB / single-crossfit vectors, which have no matrix to aggregate.

PROFILES <- list(
  base    = list(cf_variance = FALSE, tests = TRUE,  dr_rf_tests = TRUE),
  ci      = list(cf_variance = TRUE,  tests = FALSE, dr_rf_tests = FALSE),
  missing = list(cf_variance = TRUE,  tests = TRUE,  dr_rf_tests = TRUE),
  ci_mi   = list(cf_variance = TRUE,  tests = FALSE, dr_rf_tests = FALSE)
)

# ---- helpers ----------------------------------------------------------------

# grf and SuperLearner spell "weights" differently and both want NULL when absent
wts <- function(ipw, idx) if (!is.null(ipw)) ipw[idx] else NULL

# AIPW / DR pseudo-outcome
dr_pseudo <- function(Y, W, Y1.hat, Y0.hat, W.hat) {
  Y.hat <- W * Y1.hat + (1 - W) * Y0.hat
  (Y1.hat - Y0.hat) + ((Y - Y.hat) * (W - W.hat)) / (W.hat * (1 - W.hat))
}

is_binomial <- function(family) identical(family$family, "binomial")

# ---- main entry point -------------------------------------------------------

#' Fit every CATE method on one simulated dataset
#'
#' @param data data.frame laid out as the DGMs emit it: Y, W, then covariates
#' @param n_folds number of crossfitting folds V. Only the SuperLearner arm
#'   still crossfits (see the strategy note above); the others are whole-sample
#'   OOB and ignore this beyond it feeding fold_indices, which they still return.
#' @param sl_lib SuperLearner libraries, list(W = , Y = , tau = ) as
#'   R/sl_library.R's sl_libraries(n) builds them, or one character vector used
#'   for all three; NULL skips the SuperLearner arm. Learners the pretest drops
#'   are returned as `results$sl_dropped`.
#' @param fmla_info oracle formula + parameters; NULL skips the oracle arm
#' @param family gaussian() or binomial(); controls the SuperLearner outcome model
#' @param ipw optional length-n weights for the missing-data IPW arm
#' @param ci NULL, or list(boot = , sf = , alpha = ) to add half-sample bootstrap CIs
#' @param profile "base", "ci" or "missing" - see PROFILES
#' @param num.threads grf thread count, forwarded to every regression_forest()/
#'   causal_forest() call this function reaches (nuisance_rf, stage2_whole_rf,
#'   run_causal_forest, run_dr_semi_oracle). NULL (default) is grf's own default
#'   (all visible cores) - unchanged from this function's behaviour before this
#'   parameter existed. SuperLearner arms have no equivalent thread knob.
#' @param verbose_timing if TRUE, time each top-level block below with
#'   R/utils.R's `timed()` and attach the elapsed seconds as `results$timings`
#'   (a named list). Default FALSE leaves `results` exactly as before this
#'   parameter existed, so production output is unchanged.
#'
#' Each study's *_models.R is a thin shim defining `run_all_cate_methods` with
#' that study's historical signature and forwarding to this. The shared function
#' is named differently so those shims do not recurse into themselves.
#' @param Z_query optional data.frame/matrix of covariate rows (same columns,
#'   same order as data[, -c(1:2)]) to also predict CATE and, when `ci` is set,
#'   a half-sample bootstrap band at - e.g. R/dgm_scenarios.R's
#'   build_query_grid(). NULL (default) skips all of it, so every existing
#'   caller is unaffected. Adds `tau_grid` (and, with `ci`, `grid_lb`/
#'   `grid_ub`/`grid_draws`) to each arm's result list.
cate_methods <- function(data, n_folds = 10, sl_lib = NULL, fmla_info = NULL,
                         family = gaussian(), ipw = NULL, ci = NULL,
                         profile = c("base", "ci", "missing", "ci_mi"),
                         num.threads = NULL, verbose_timing = FALSE,
                         Z_query = NULL) {

  profile <- match.arg(profile)
  p <- PROFILES[[profile]]
  sl_lib <- as_sl_libs(sl_lib)

  X <- as.matrix(data[, -c(1:2)])
  Y <- data$Y
  W <- data$W
  n_obs <- nrow(X)

  Z_query_mat <- if (!is.null(Z_query)) as.matrix(Z_query) else NULL

  fold_indices <- sort(seq(n_obs) %% n_folds) + 1
  fold_list <- unique(fold_indices)

  results <- list()
  timings <- list()

  # expr is an R promise, evaluated exactly once on first access whichever
  # branch runs below - so this never changes what gets computed, only whether
  # the elapsed time is captured alongside it.
  time_step <- function(name, expr) {
    if (verbose_timing) {
      t <- timed(expr)
      timings[[name]] <<- t$time
      t$value
    } else {
      expr
    }
  }

  cat("Computing nuisance functions...\n")
  nuisances_rf <- time_step("nuisance_rf", nuisance_rf(X, Y, W, ipw, num.threads = num.threads))

  cat("Running Causal Forest...\n")
  results$causal_forest <- time_step("causal_forest",
    run_causal_forest(X, Y, W, nuisances_rf, ipw,
                      variance = p$cf_variance,
                      tests = p$tests, num.threads = num.threads, Z_query = Z_query_mat))
  if (!is.null(ci)) {
    cat("Running Causal Forest bootstrap... \n")
    results$causal_forest <- c(results$causal_forest,
      cf_oob_half_boot(X, Y, W, results$causal_forest, results$causal_forest$tau,
                       ci$boot, ci$sf, ci$alpha,
                       Z_query = Z_query_mat, tau_grid = results$causal_forest$tau_grid))
  }

  cat("Running DR Random Forest...\n")
  results$dr_random_forest <- time_step("dr_random_forest", if (p$dr_rf_tests) {
    run_dr_random_forest(X, Y, W, nuisances_rf, ipw, num.threads = num.threads, Z_query = Z_query_mat)
  } else {
    s <- stage2_whole_rf(X, nuisances_rf$po, ipw, num.threads = num.threads, Z_query = Z_query_mat)
    list(tau = s$tau, tau_grid = s$tau_grid)
  })
  if (!is.null(ci)) {
    cat("Running DR RF bootstrap... \n")
    results$dr_random_forest <- c(results$dr_random_forest,
      rf_oob_half_boot(X, Y, W, nuisances_rf$po, results$dr_random_forest$tau,
                       ci$boot, ci$sf, ci$alpha,
                       Z_query = Z_query_mat, tau_grid = results$dr_random_forest$tau_grid))
  }

  # the missing-data variant skips the arms that need a complete covariate matrix
  complete_X <- !anyNA(X)

  if (!is.null(fmla_info) && (profile != "missing" || complete_X)) {
    cat("Running DR Oracle...\n")
    results$dr_oracle <- time_step("dr_oracle",
      run_dr_oracle(X, Y, W, fmla_info, ipw,
                    tests = p$tests, num.threads = num.threads,
                    Z_query = Z_query_mat))
    if (!is.null(ci)) {
      cat("Runnings Oracle bootstrap...\n")
      results$dr_oracle <- c(results$dr_oracle,
        rf_oob_half_boot(X, Y, W, results$dr_oracle$po, results$dr_oracle$tau,
                         ci$boot, ci$sf, ci$alpha,
                         Z_query = Z_query_mat, tau_grid = results$dr_oracle$tau_grid))
    }
  }

  if (profile != "missing" || complete_X) {
    cat("Running DR Semi-Oracle...\n")
    results$dr_semi_oracle <- time_step("dr_semi_oracle",
      run_dr_semi_oracle(X, Y, W, ipw, tests = p$tests, num.threads = num.threads,
                        Z_query = Z_query_mat))
    if (!is.null(ci)) {
      cat("Running Semi-Oracle bootstrap...\n")
      results$dr_semi_oracle <- c(results$dr_semi_oracle,
        rf_oob_half_boot(X, Y, W, results$dr_semi_oracle$po, results$dr_semi_oracle$tau,
                         ci$boot, ci$sf, ci$alpha,
                         Z_query = Z_query_mat, tau_grid = results$dr_semi_oracle$tau_grid))
    }
  }

  if (!is.null(sl_lib) && (profile != "missing" || complete_X)) {
    cat("Running DR SuperLearner...\n")
    X <- as.data.frame(X)
    results$dr_superlearner <- time_step("dr_superlearner", {
      nuisances_sl <- nuisance_sl(X, Y, W, fold_indices, sl_lib, ipw, family = family)
      out <- run_dr_superlearner(X, Y, W, nuisances_sl,
                                 fold_indices, fold_list,
                                 sl_lib, ipw, tests = p$tests)
      # this block is a promise forced inside time_step(), but it was created
      # (and so evaluates) in cate_methods' own frame - plain `<-` reaches
      # cate_methods' `results` directly; `<<-` here would skip that frame and
      # write into whatever encloses cate_methods instead.
      # The pretest's dropped learners ride along as attributes; they are
      # saved on their own so the nuisance and tau fields keep their shape.
      results$sl_dropped <- rbind(attr(nuisances_sl, "sl_dropped"), out$sl_dropped)
      attr(nuisances_sl, "sl_dropped") <- NULL
      out$sl_dropped <- NULL
      results$nuisances_sl <- nuisances_sl
      out
    })
  }

  results$nuisances_rf <- nuisances_rf
  results$fold_indices <- fold_indices
  if (verbose_timing) results$timings <- timings

  results
}

# ---- stage 1: nuisance estimation -------------------------------------------

#' One arm's outcome forest: a regression forest on the rows with W == arm
#'
#' The per-arm fit shared by t_learner_rf (whole sample) and
#' t_learner_rf_split (one train/test split), so the two differ only in which
#' rows they are handed and how they predict.
arm_forest <- function(X, Y, W, arm, ipw = NULL, num.threads = NULL) {
  in_arm <- W == arm
  regression_forest(X[in_arm, , drop = FALSE], Y[in_arm],
                    sample.weights = wts(ipw, in_arm),
                    num.threads = num.threads)
}

#' T-learner outcome model with regression forests: one forest per arm
#'
#' Whole-sample, no splitting. A unit's prediction from its own arm's forest is
#' out-of-bag; its prediction from the other arm's forest is an ordinary
#' newdata prediction from a forest that never saw it - so both are honest.
#' This is the "oob_oob" arm of crossfitting/cf_models.R, which calls it
#' (through nuisance_rf) rather than keeping its own copy.
#' @return list(Y0.hat, Y1.hat), each length n
t_learner_rf <- function(X, Y, W, ipw = NULL, num.threads = NULL) {
  n_obs <- nrow(X)
  arm_fit <- function(arm) {
    in_arm <- W == arm
    forest <- arm_forest(X, Y, W, arm, ipw, num.threads)
    pred <- numeric(n_obs)
    pred[in_arm] <- predict(forest)$predictions
    pred[!in_arm] <- predict(forest, newdata = X[!in_arm, , drop = FALSE])$predictions
    pred
  }
  # control arm first: fixes the order the two fits consume the RNG stream
  Y0.hat <- arm_fit(0)
  Y1.hat <- arm_fit(1)
  list(Y0.hat = Y0.hat, Y1.hat = Y1.hat)
}

#' T-learner outcome model with regression forests, for one train/test split
#'
#' One forest per arm, each fit on that arm's training rows only and predicting
#' every test row. No production estimator crossfits its forests any more; this
#' is for crossfitting/cf_models.R's crossfit arms (dcf, scf_*), which then
#' differ from the whole-sample t_learner_rf only in splitting.
#' @param in_train,in_test logical row masks
#' @return list(Y0.hat, Y1.hat), each of length sum(in_test)
t_learner_rf_split <- function(X, Y, W, in_train, in_test, ipw = NULL,
                               num.threads = NULL) {
  X_train <- X[in_train, , drop = FALSE]
  X_test <- X[in_test, , drop = FALSE]
  ipw_train <- if (!is.null(ipw)) ipw[in_train] else NULL
  arm_pred <- function(arm) {
    forest <- arm_forest(X_train, Y[in_train], W[in_train], arm, ipw_train, num.threads)
    predict(forest, newdata = X_test)$predictions
  }
  # control arm first, as in t_learner_rf
  Y0.hat <- arm_pred(0)
  Y1.hat <- arm_pred(1)
  list(Y0.hat = Y0.hat, Y1.hat = Y1.hat)
}

#' Whole-sample OOB nuisance estimation with regression forests (T-learner)
#'
#' No sample splitting: t_learner_rf supplies Y0.hat/Y1.hat, and two more
#' whole-sample forests supply W.hat and Y.hat.cf, both taken out-of-bag. This
#' is the "oob_oob" arm of crossfitting/cf_models.R. It replaced the S-learner
#' "oob_oob_s" arm (one forest on cbind(W, X), read at counterfactual rows via
#' grf's X.orig) when the DR-learners moved to per-arm outcome models; that
#' arm had itself replaced the double-crossfit this function used to do.
nuisance_rf <- function(X, Y, W, ipw = NULL, num.threads = NULL) {

  mu <- t_learner_rf(X, Y, W, ipw, num.threads = num.threads)
  Y0.hat <- mu$Y0.hat
  Y1.hat <- mu$Y1.hat

  W.hat <- trim_ps(predict(regression_forest(X, W, sample.weights = ipw,
                                             num.threads = num.threads))$predictions)
  Y.hat.cf <- predict(regression_forest(X, Y, sample.weights = ipw,
                                        num.threads = num.threads))$predictions

  Y.hat <- W * Y1.hat + (1 - W) * Y0.hat
  po <- dr_pseudo(Y, W, Y1.hat, Y0.hat, W.hat)

  list(po = po, Y.hat = Y.hat, Y.hat.cf = Y.hat.cf, Y0.hat = Y0.hat, W.hat = W.hat)
}

#' One train/test split's SuperLearner nuisances (T-learner outcome model)
#'
#' One SuperLearner per arm on X, each fit on that arm's training rows only
#' and predicting every test row, plus a propensity SuperLearner on all the
#' training rows. The per-split body of nuisance_sl, shared with
#' crossfitting/cf_models.R's double-crossfit arm so that its fold-pair fits
#' are the production estimator's.
#'
#' @param X covariates, as a data frame
#' @param in_train,in_test logical row masks
#' @param sl_lib list(W = , Y = , tau = ), already through as_sl_libs()
#' @return list(po, Y.hat, Y0.hat, W.hat) at the test rows, and libs - the
#'   pretested libraries list(Y0 = , Y1 = , W = ) for dropped_table()
sl_split_fit <- function(X, Y, W, in_train, in_test, sl_lib, ipw = NULL,
                         family = gaussian()) {

  binom <- is_binomial(family)

  X_train <- X[in_train, ]
  X_test <- X[in_test, ]
  Y_family <- if (binom) binomial() else gaussian()

  # one outcome model per arm, each predicting at every held-out row. The
  # failsafe - SuperLearner returns all-zero predictions when every learner
  # ends up with zero weight - falls back to that arm's training mean. A fit
  # that errors outright falls back to the mean inside sl_fit_predict, and
  # is recorded with the pretest's drops.
  arm_fit <- function(arm) {
    in_arm <- in_train & W == arm
    lib <- pretest_superlearner(Y[in_arm], X[in_arm, ], sl_lib$Y, Y_family)
    fit <- sl_fit_predict(Y[in_arm], X[in_arm, ], list(test = X_test), lib,
                          family = Y_family, obsWeights = wts(ipw, in_arm))
    pred <- fit$pred$test
    if (all(pred == 0)) {
      warning("SuperLearner failed for Y.hat in arm W = ", arm, ". Using its mean.")
      pred <- rep(mean(Y[in_arm], na.rm = TRUE), sum(in_test))
    }
    list(pred = pred, lib = mark_failed_fit(lib, fit))
  }
  fit0 <- arm_fit(0)
  fit1 <- arm_fit(1)

  W_lib <- pretest_superlearner(W[in_train], X_train, sl_lib$W, binomial())
  W_fit <- sl_fit_predict(W[in_train], X_train, list(w = X_test), W_lib,
                          family = binomial(), obsWeights = wts(ipw, in_train))
  W_lib <- mark_failed_fit(W_lib, W_fit)
  W.hat <- W_fit$pred$w
  if (all(W.hat == 0)) {
    warning("SuperLearner failed for W.hat. Using mean(W).")
    W.hat <- rep(mean(W[in_train], na.rm = TRUE), sum(in_test))
  }

  # clamp propensities away from 0/1, same as nuisance_rf's W.hat
  W.hat <- trim_ps(W.hat)

  Y0.hat <- fit0$pred
  Y1.hat <- fit1$pred
  W_test <- W[in_test]
  Y.hat <- W_test * Y1.hat + (1 - W_test) * Y0.hat
  po <- dr_pseudo(Y[in_test], W_test, Y1.hat, Y0.hat, W.hat)

  list(po = po, Y.hat = Y.hat, Y0.hat = Y0.hat, W.hat = W.hat,
       libs = list(Y0 = fit0$lib, Y1 = fit1$lib, W = W_lib))
}

#' Single leave-one-fold-out nuisance estimation with SuperLearner
#'
#' One split, shared with the stage-2 regression via the same fold_indices
#' (see run_dr_superlearner / stage_2_sl) rather than double-crossfit
#' nuisances feeding a separately-split stage 2 - the "scf_scf" arm of
#' crossfitting/cf_models.R, which calls this directly. The outcome model is a
#' T-learner - see sl_split_fit.
#'
#' @param sl_lib list(W = , Y = , tau = ) - see as_sl_libs(). The pretest's
#'   dropped learners are returned as attr(, "sl_dropped").
nuisance_sl <- function(X, Y, W, fold_indices, sl_lib, ipw = NULL,
                        family = gaussian()) {

  sl_lib <- as_sl_libs(sl_lib)

  cross_fits <- future_map(unique(fold_indices), function(fold) {
    in_train <- fold_indices != fold
    fit <- sl_split_fit(X, Y, W, in_train, !in_train, sl_lib, ipw, family)
    c(list(fold = fold), fit[c("po", "Y.hat", "Y0.hat", "W.hat")],
      list(dropped = dropped_table(fit$libs, fold)))
  }, .options = furrr_options(seed = TRUE))

  out <- scatter_folds(cross_fits, fold_indices, c("po", "Y.hat", "Y0.hat", "W.hat"))
  attr(out, "sl_dropped") <- do.call(rbind, lapply(cross_fits, `[[`, "dropped"))
  out
}

# pretest_superlearner() lives in R/sl_library.R with the libraries it tests.

# ---- stage 2: final CATE regression -----------------------------------------

#' Whole-sample OOB second stage: one forest, its OOB predictions
#'
#' var_oob is grf's own OOB variance estimate (bootstrap of little bags), free
#' alongside the predictions since regression_forest already defaults to
#' ci.group.size = 2. predict() does not consume R's RNG stream, so asking for
#' the variance leaves the point estimate unchanged. Ported from
#' crossfitting/cf_models.R::stage2_whole_rf, replacing the leave-one-fold-out
#' stage_2_rf this function used to be for dr_random_forest, dr_oracle and
#' dr_semi_oracle alike.
#' @param Z_query optional covariate rows (matrix/data.frame, same columns as
#'   X) to also predict at, off this same fitted forest, before it goes out of
#'   scope. NULL (default) adds nothing - see cate_methods' Z_query doc.
stage2_whole_rf <- function(X, po, ipw = NULL, num.threads = NULL, Z_query = NULL) {
  forest <- regression_forest(X, po, sample.weights = ipw, num.threads = num.threads)
  pred <- predict(forest, estimate.variance = TRUE)
  tau_grid <- if (!is.null(Z_query)) predict(forest, newdata = Z_query)$predictions else NULL
  list(tau = pred$predictions, variance = pred$variance.estimates, tau_grid = tau_grid)
}

#' Crossfit second stage with SuperLearner
#'
#' @param sl_lib list(W = , Y = , tau = ) or one character vector - see
#'   as_sl_libs(); the tau library is used. The pretest's dropped learners are
#'   returned as attr(tau, "sl_dropped").
stage_2_sl <- function(X, po, fold_indices, fold_list, sl_lib, ipw = NULL) {
  n_obs <- nrow(X)
  single <- is.vector(po)
  tau_lib <- as_sl_libs(sl_lib)$tau

  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold <- !in_train

    y_train <- if (single) po[in_train] else po[in_train, fold]
    po_lib <- pretest_superlearner(y_train, X[in_train, ], tau_lib, gaussian())
    po_fit <- sl_fit_predict(y_train, X[in_train, ], list(tau = X[in_fold, ]), po_lib,
                             family = gaussian(), obsWeights = wts(ipw, in_train))
    list(fold = fold, predictions = po_fit$pred$tau,
         dropped = dropped_table(list(tau = mark_failed_fit(po_lib, po_fit)), fold))
  }, .options = furrr_options(seed = TRUE))

  tau <- rep(NA, n_obs)
  for (result in tau_results) tau[fold_indices == result$fold] <- result$predictions
  attr(tau, "sl_dropped") <- do.call(rbind, lapply(tau_results, `[[`, "dropped"))
  tau
}

# ---- the estimators ---------------------------------------------------------

#' Causal forest using grf's own internal cross-fitting
#'
#' No externally-supplied nuisances: leaving Y.hat/W.hat NULL makes grf
#' cross-fit them internally, and tau is grf's own out-of-bag prediction. This
#' is the "cf_default" arm validated in crossfitting/cf_models.R against the
#' fold-wise external-crossfit alternative this function used to implement.
#'
#' `nuisances` is only used for the BLP/independence tests below - it is the
#' nuisance_rf() object shared with dr_random_forest, not what fits the forest.
#' The forest's own Y.hat/W.hat are returned (as Y.hat.cf/W.hat, matching the
#' field-naming convention) so the half-sample bootstrap can hold them fixed.
run_causal_forest <- function(X, Y, W, nuisances, ipw = NULL, variance = FALSE,
                              tests = TRUE, num.threads = NULL, Z_query = NULL) {
  forest <- causal_forest(X, Y, W, sample.weights = ipw, num.threads = num.threads)
  pred <- predict(forest, estimate.variance = variance)

  out <- list(tau = pred$predictions, Y.hat.cf = forest$Y.hat, W.hat = forest$W.hat)
  if (variance) out$variance <- pred$variance.estimates
  if (!is.null(Z_query)) out$tau_grid <- predict(forest, newdata = Z_query)$predictions
  if (tests) {
    out$BLP_whole <- run_blp_whole(Y, W, nuisances$W.hat, nuisances$Y0.hat, out$tau)
    out$independence_cate <- run_independence_test_whole(X, out$tau)
    out$independence_po <- run_independence_test_whole(X, nuisances$po)
  }
  out
}

#' DR-learner with a whole-sample OOB regression-forest second stage
run_dr_random_forest <- function(X, Y, W, nuisances, ipw = NULL, tests = TRUE,
                                 num.threads = NULL, Z_query = NULL) {
  s <- stage2_whole_rf(X, nuisances$po, ipw, num.threads = num.threads, Z_query = Z_query)

  out <- list(tau = s$tau, variance = s$variance, tau_grid = s$tau_grid)
  if (tests) {
    out$BLP_whole <- run_blp_whole(Y, W, nuisances$W.hat, nuisances$Y0.hat, out$tau)
    out$independence_cate <- run_independence_test_whole(X, out$tau)
    out$independence_po <- run_independence_test_whole(X, nuisances$po)
  }
  out
}

#' DR-learner with the true outcome model, a known propensity of 0.5, and a
#' whole-sample OOB second stage
#'
#' fmla_info$fmla returns the outcome mean E[Y | X, W] itself - see the header.
run_dr_oracle <- function(X, Y, W, fmla_info, ipw = NULL, tests = TRUE,
                          num.threads = NULL, Z_query = NULL) {
  n_obs <- nrow(X)

  X <- as.data.frame(X)
  list2env(fmla_info$params, envir = environment())
  fmla <- parse(text = fmla_info$fmla)

  W_temp <- rep(1, n_obs)
  Y1.hat <- eval(fmla, envir = list2env(c(list(W = W_temp), X)))

  W_temp <- rep(0, n_obs)
  Y0.hat <- eval(fmla, envir = list2env(c(list(W = W_temp), X)))

  Y.hat <- eval(fmla, envir = list2env(c(list(W = W), X)))
  W.hat <- rep(0.5, n_obs)

  X <- as.matrix(X)

  po <- (Y1.hat - Y0.hat) + ((Y - Y.hat) * (W - W.hat)) / (W.hat * (1 - W.hat))
  s <- stage2_whole_rf(X, po, ipw, num.threads = num.threads, Z_query = Z_query)

  out <- list(tau = s$tau, variance = s$variance, tau_grid = s$tau_grid, po = po, Y0.hat = Y0.hat)
  if (tests) {
    out$BLP_whole <- run_blp_whole(Y, W, W.hat, Y0.hat, out$tau)
    out$independence_cate <- run_independence_test_whole(X, out$tau)
    out$independence_po <- run_independence_test_whole(X, po)
  }
  out
}

#' DR-learner with a known propensity of 0.5, a whole-sample OOB outcome
#' model, and a whole-sample OOB second stage
run_dr_semi_oracle <- function(X, Y, W, ipw = NULL, tests = TRUE, num.threads = NULL,
                               Z_query = NULL) {
  n_obs <- nrow(X)
  W.hat <- rep(0.5, n_obs)

  # the same per-arm outcome forests as nuisance_rf, so the gap to
  # dr_random_forest is the propensity alone
  mu <- t_learner_rf(X, Y, W, ipw, num.threads = num.threads)
  Y0.hat <- mu$Y0.hat
  Y1.hat <- mu$Y1.hat

  po <- dr_pseudo(Y, W, Y1.hat, Y0.hat, W.hat)
  s <- stage2_whole_rf(X, po, ipw, num.threads = num.threads, Z_query = Z_query)

  out <- list(tau = s$tau, variance = s$variance, tau_grid = s$tau_grid, po = po, Y0.hat = Y0.hat)
  if (tests) {
    out$BLP_whole <- run_blp_whole(Y, W, W.hat, Y0.hat, out$tau)
    out$independence_cate <- run_independence_test_whole(X, out$tau)
    out$independence_po <- run_independence_test_whole(X, po)
  }
  out
}

#' DR-learner with SuperLearner nuisances and second stage, sharing one split
run_dr_superlearner <- function(X, Y, W, nuisances, fold_indices, fold_list,
                                sl_lib, ipw = NULL, tests = TRUE) {
  tau <- stage_2_sl(X, nuisances$po, fold_indices, fold_list, sl_lib, ipw)
  sl_dropped <- attr(tau, "sl_dropped")
  attr(tau, "sl_dropped") <- NULL

  out <- list(tau = tau, sl_dropped = sl_dropped)
  if (tests) {
    out$BLP_whole <- run_blp_whole(Y, W, nuisances$W.hat, nuisances$Y0.hat, tau)
    out$independence_cate <- run_independence_test_whole(X, tau)
    out$independence_po <- run_independence_test_whole(X, nuisances$po)
  }
  out
}

# ---- post-estimation heterogeneity tests ------------------------------------

#' Is x constant, up to floating-point rounding?
#'
#' Both tests below are meaningless for a constant CATE, and neither fails
#' cleanly on one. coin::independence_test() only warns ("zero diagonal
#' elements") and returns a huge statistic with p ~ 0 - a spurious rejection.
#' GenericML::BLP() errors only when tau is EXACTLY constant (bug L below). A
#' true CATE can be constant up to rounding instead: truth_at() computes tau as
#' (p0 + bW) - p0, which leaves noise of ~1e-17 in scenario 1, and BLP then fits
#' beta.2 ~ 1e15 on it. The tolerance is relative, far below any real variation
#' in an estimated or true CATE.
is_constant <- function(x, tol = sqrt(.Machine$double.eps)) {
  x <- x[!is.na(x)]
  length(x) < 2 || sd(x) <= tol * max(1, mean(abs(x)))
}

# Best Linear Predictor of the CATE (GenericML). Returns the whole coefficient
# block - Estimate, Std. Error, t value, Pr(>|t|) - with the residual df as
# attr(, "df"); blp_p_value() below reads beta.2's p-value from it. It used to
# keep only columns 1 and 4, which dropped the standard error a
# multiple-imputation pooling rule needs (mi_test_table() below).
# blp_p_value() reads both shapes, so results saved before the change still
# read correctly.
#
# vcov_type is the sandwich::vcovHC type. The default "const" (homoskedastic
# OLS SEs, GenericML's own default) is what every estimation-time call uses, so
# the saved BLP_whole objects are unchanged; R/metrics.R::hte_test_metrics()
# recomputes with "HC3" for BLP_p_os - see blp_inputs().
#
# bug L: GenericML::BLP() regresses on beta.2 = (W - W.hat) * (tau - mean(tau)).
# When tau is exactly constant (a degenerate/near-constant CATE fit - seen with
# scenario 4's low-amplitude cos(X4) effect at small n, especially once
# pretest_superlearner has whittled a fold's library down to 1-2 survivors),
# beta.2 becomes identically zero, lm() marks it aliased (coef = NA), and
# sandwich::vcovHC() drops that coefficient's row/column entirely rather than
# keeping it as NA - so GenericML's internal indexing by name throws "subscript
# out of bounds". NULL is a deliberate fallback, not just a safe default:
# hte_test_metrics() (R/metrics.R) already maps BLP_whole = NULL to BLP_p = NA,
# and NA is the statistically correct answer here - beta.2 has no fitted
# coefficient to attach a p-value to when tau has zero variance. is_constant()
# extends the same fallback to a tau that is constant only up to rounding.
run_blp_whole <- function(Y, W, W.hat, Y0.hat, tau, vcov_type = "const") {
  if (is_constant(tau)) {
    warning("run_blp_whole: tau is constant (up to rounding); returning NULL.")
    return(NULL)
  }
  vcov_control <- setup_vcov(estimator = "vcovHC", arguments = list(type = vcov_type))
  tryCatch(
    # unclass: a plain matrix, so reading it back needs no lmtest; the df
    # attribute survives
    unclass(BLP(Y, W, W.hat, Y0.hat, tau, vcov_control = vcov_control)$coefficients),
    error = function(e) {
      warning("run_blp_whole: BLP() failed (likely a constant/degenerate tau); ",
              "returning NULL. ", conditionMessage(e))
      NULL
    }
  )
}

#' beta.2's p-value from a run_blp_whole() coefficient block
#'
#' "two" is Pr(>|t|), what BLP_p has always been. "one" tests H1: beta.2 > 0,
#' the direction Chernozhukov et al. and grf::test_calibration() use: a proxy
#' that ranks units in reverse (beta.2 < 0) is not evidence of heterogeneity,
#' but it gets the same two-sided p as a correct one.
#'
#' Reads both block shapes run_blp_whole() has saved: the full block
#' (Estimate, Std. Error, t value, Pr(>|t|), attr "df") and the older
#' Estimate + Pr(>|t|) pair, for which the one-sided p is recovered from the
#' two-sided one and the sign of the estimate.
#' @param blp a run_blp_whole() result; NULL (degenerate tau) gives NA
blp_p_value <- function(blp, sided = c("one", "two")) {
  sided <- match.arg(sided)
  if (is.null(blp)) return(NA_real_)
  p_two <- unname(blp["beta.2", ncol(blp)])
  if (sided == "two") return(p_two)
  est <- unname(blp["beta.2", 1])
  if (ncol(blp) < 4) {
    return(if (est > 0) p_two / 2 else 1 - p_two / 2)
  }
  stat <- est / unname(blp["beta.2", 2])
  df <- attr(blp, "df")
  if (is.null(df)) pnorm(stat, lower.tail = FALSE) else pt(stat, df, lower.tail = FALSE)
}

#' The inputs one model's BLP test was run on, rebuilt from a saved run
#'
#' Lets R/metrics.R::hte_test_metrics() recompute the BLP (with HC3 SEs, for
#' BLP_p_os) from results already on disk. MUST stay in step with the
#' run_blp_whole() calls in the run_* functions above:
#'   causal_forest, dr_random_forest  nuisances_rf's W.hat and Y0.hat
#'   dr_oracle, dr_semi_oracle        W.hat = 0.5, the arm's own Y0.hat
#'   dr_superlearner                  nuisances_sl's W.hat and Y0.hat
#' and with R/metrics.R::add_t_learners(), whose T-learners reuse their DR
#' counterpart's:
#'   t_random_forest                  nuisances_rf's W.hat and Y0.hat
#'   t_superlearner                   nuisances_sl's W.hat and Y0.hat
#' @param sim_res one run's saved results object
#' @param model model name
#' @return list(Y, W, W.hat, Y0.hat, tau), or NULL when the run has no single
#'   dataset (multiple_imputation saves a list of them) or lacks a field
blp_inputs <- function(sim_res, model) {
  if (!is.data.frame(sim_res$data)) return(NULL)
  m <- sim_res[[model]]
  nuis <- switch(model,
    causal_forest = , dr_random_forest = , t_random_forest = sim_res$nuisances_rf,
    dr_oracle = , dr_semi_oracle = list(W.hat = 0.5, Y0.hat = m$Y0.hat),
    dr_superlearner = , t_superlearner = sim_res$nuisances_sl,
    NULL
  )
  if (is.null(m$tau) || is.null(nuis$W.hat) || is.null(nuis$Y0.hat)) return(NULL)
  n_obs <- nrow(sim_res$data)
  list(Y = sim_res$data$Y, W = sim_res$data$W,
       W.hat = rep_len(nuis$W.hat, n_obs), Y0.hat = nuis$Y0.hat, tau = m$tau)
}

# Omnibus independence test of the estimated CATEs against the covariates.
# The test is the one WATCH (github.com/Novartis/WATCH) runs on the DR
# pseudo-outcome, and is kept exactly as is. A constant tau returns NA before
# coin is called: coin does not fail on one, it only warns and returns p ~ 0
# (see is_constant()).
run_independence_test_whole <- function(X, tau) {
  if (is_constant(tau)) {
    return(list(p_value = NA_real_, statistic = NA_real_, df = NA_real_,
                method = "constant_input"))
  }
  test_data <- data.frame(tau = tau, X)
  tryCatch({
    test_result <- coin::independence_test(
      tau ~ .,
      data = test_data,
      teststat = "quadratic"
    )
    list(
      p_value = coin::pvalue(test_result),
      statistic = coin::statistic(test_result),
      # the asymptotic chi-square df - what a pooling rule on the statistic
      # (e.g. D2) needs alongside it; see mi_test_table()
      df = test_result@statistic@df,
      method = "independence_test"
    )
  }, error = function(e) {
    list(p_value = 1, statistic = 0, df = NA_real_,
         method = "independence_test_failed")
  })
}

#' BLP and independence tests run on the true CATE and true nuisances
#'
#' The oracle counterpart of the per-estimator BLP_whole / independence_cate
#' fields: same test functions, same call shape, but tau = true CATE,
#' Y0.hat = true p0, W.hat = 0.5 (known exactly - see generate_scenario_data(),
#' which draws W <- rbinom(n, 1, 0.5) unconditionally). No independence_po
#' counterpart: dr_oracle's own independence_po already tests the true
#' pseudo-outcome against X.
#'
#' Scenario 1's true CATE is constant up to rounding, so both tests come back
#' NA there (is_constant()) - the true-CATE rows have no null scenario.
#' BLP_whole_hc3 is the same BLP with HC3 SEs, for BLP_p_os.
run_true_cate_tests <- function(X, Y, W, truth) {
  W.hat <- rep(0.5, length(W))
  list(
    BLP_whole = run_blp_whole(Y, W, W.hat, truth$p0, truth$tau),
    BLP_whole_hc3 = run_blp_whole(Y, W, W.hat, truth$p0, truth$tau, vcov_type = "HC3"),
    independence_cate = run_independence_test_whole(X, truth$tau)
  )
}

#' One row of true-CATE HTE test p-values for one run
#'
#' NA row when sim_res$data is not a single data.frame - true today only for
#' multiple_imputation rows (missing/binary, missing/continuous), which save
#' `data` as a list of 50 imputed data.frames with no single X to test
#' against. Same open pooling question as the estimated-CATE tests - see
#' mi_test_table() and missing/README.md.
true_cate_test_row <- function(sim_res) {
  if (!is.data.frame(sim_res$data)) {
    return(tibble::tibble(BLP_p = NA_real_, BLP_p_os = NA_real_,
                          indep_cate = NA_real_))
  }
  X <- as.matrix(sim_res$data[, -c(1:2)])
  Y <- sim_res$data$Y
  W <- sim_res$data$W
  out <- run_true_cate_tests(X, Y, W, sim_res$truth)
  tibble::tibble(
    BLP_p = blp_p_value(out$BLP_whole, "two"),
    BLP_p_os = blp_p_value(out$BLP_whole_hc3, "one"),
    indep_cate = as.numeric(out$independence_cate$p_value)
  )
}

# ---- multiple imputation ----------------------------------------------------

#' Rubin-combine CATE estimates across multiply-imputed datasets
#'
#' @param res_list one run_all_cate_methods result per imputation
#' @param model which model's estimates to combine
combine_mi <- function(res_list, model) {
  res <- list()

  tau_list <- lapply(res_list, function(x) x[[model]][["tau"]])
  tau_mat <- do.call(cbind, tau_list)
  res$tau <- rowMeans(tau_mat)

  var_list <- lapply(res_list, function(x) x[[model]][["variance"]])
  var_mat <- do.call(cbind, var_list)
  if (!is.null(var_mat)) {
    w_var <- rowMeans(var_mat)
    b_var <- apply(tau_mat, 1, var)
    res$variance <- w_var + (1 + 1 / length(res_list)) * b_var
  }

  res
}

#' Per-imputation HTE test results for one model, unpooled
#'
#' Every imputation's fit already runs the BLP and independence tests (profile
#' "missing" has tests on), but combine_mi() pools only tau and variance, so
#' they used to be discarded - which is why the multiple_imputation arm had no
#' HTE tests. This keeps what any candidate pooling rule needs, one row per
#' imputation: beta.2's estimate, SE and residual df (Rubin's rules), each
#' independence test's chi-square statistic and df (D2), and every p-value (a
#' p-value combination rule).
#'
#' No rule is applied. Which one is still an open decision (missing/README.md),
#' and saving the rows means it can be made at metrics time without a re-run.
#' Until then the arm has no BLP_whole / independence_* fields, so
#' hte_test_metrics() still reports NA for it.
#'
#' @param res_list one run_all_cate_methods result per imputation
#' @param model which model's tests to extract
mi_test_table <- function(res_list, model) {
  bind_rows(lapply(seq_along(res_list), function(k) {
    m <- res_list[[k]][[model]]
    blp <- m$BLP_whole
    # NULL is bug L's degenerate-tau fallback (run_blp_whole) - kept as an NA
    # row rather than dropped, so the table always has one row per imputation
    # by position, not name: a coeftest block is always Estimate, Std. Error,
    # statistic, p-value, but the last two are named t or z depending on the
    # vcov, and a name lookup failing here would lose the whole 50-imputation run
    blp_col <- function(col) if (is.null(blp)) NA_real_ else unname(blp["beta.2", col])
    # a failed independence test comes back as p = 1, statistic = 0, and a
    # constant tau as p = NA (see run_independence_test_whole); df = NA marks both
    indep <- function(test, field) {
      v <- test[[field]]
      if (is.null(v)) NA_real_ else as.numeric(v)
    }
    tibble::tibble(
      imputation      = k,
      blp_estimate    = blp_col(1),
      blp_se          = blp_col(2),
      blp_df          = if (is.null(blp) || is.null(attr(blp, "df"))) {
        NA_real_
      } else as.numeric(attr(blp, "df")),
      blp_p           = blp_col(4),
      indep_cate_stat = indep(m$independence_cate, "statistic"),
      indep_cate_df   = indep(m$independence_cate, "df"),
      indep_cate_p    = indep(m$independence_cate, "p_value"),
      indep_po_stat   = indep(m$independence_po, "statistic"),
      indep_po_df     = indep(m$independence_po, "df"),
      indep_po_p      = indep(m$independence_po, "p_value")
    )
  }))
}
