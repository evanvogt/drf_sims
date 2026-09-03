##########
# title: shared CATE estimators (DR-learner family + causal forest)
##########
# One implementation of what used to be seven near-identical copies:
#   continuous/cts_models.R              binary/bin_models.R
#   missing/continuous/cts_miss_models.R missing/binary/bin_miss_models.R
#   missing/ci_example/cts_miss_ci_models.R
#   confidence_intervals/continuous/cts_ci_models.R
#   confidence_intervals/binary/bin_ci_models.R
#
# Those copies differed on four axes, which are the four arguments below:
#
#   family   gaussian() vs binomial(), controlling the SuperLearner outcome model
#            (family + method.NNloglik).
#   oracle_link  where the inverse link lives for the oracle arm. The repo has
#            TWO conventions and they must not be conflated with `family`:
#              binary/, confidence_intervals/binary/  -> get_binary_oracle_info
#                returns a LINEAR PREDICTOR formula and the model code applies
#                plogis. oracle_link = "logit".
#              missing/binary/  -> bin_miss_dgms.R bakes plogis(...) into the
#                formula string itself and the model applies none.
#                oracle_link = "identity".
#            Both give the same answer; applying plogis to a formula that already
#            contains it does not. Worth harmonising when the DGMs are unified.
#   ipw      sample.weights (grf) / obsWeights (SuperLearner) for the missing-data
#            IPW arm. NULL reproduces the unweighted path exactly.
#   ci       list(boot=, sf=, alpha=) turns on the half-sample bootstrap.
#   profile  which historical variant's *orchestration* to reproduce - see below.
#
# `profile` exists because the variants also disagree about which post-estimation
# tests get run, and those disagreements look like drift rather than design. They
# are reproduced exactly rather than harmonised, because harmonising would change
# results for studies that are not otherwise re-running. See PROFILES.
#
# Crossfitting strategy, per crossfitting/'s comparison of alternatives against
# the double-crossfitting this file used to do throughout:
#   dr_random_forest, dr_oracle, dr_semi_oracle  whole-sample OOB, "oob_oob_s":
#     an S-learner nuisance forest (nuisance_rf) with no sample splitting, OOB
#     counterfactual predictions via oob_predict_counterfactual, and an OOB
#     stage-2 regression forest (stage2_whole_rf).
#   causal_forest   grf's own internal cross-fitting, "cf_default": a plain
#     causal_forest(X, Y, W) with no externally-supplied nuisances.
#   dr_superlearner   single leave-one-fold-out crossfit, "scf_scf": nuisance_sl
#     and stage_2_sl share the same fold_indices, rather than double-crossfit
#     nuisances feeding a separately-split stage 2.
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
# The existing ../results/missing/{binary,continuous} files were NOT re-run for
# it. Both tests are deterministic (GenericML::BLP is an OLS with a sandwich
# vcov; coin::independence_test is asymptotic under teststat = "quadratic"), and
# the saved results retain nuisances_rf, data and tau, so the three fields were
# recomputed exactly in place by R/patch_hte_tests.R. One asymmetry follows and
# is deliberate: run_dr_random_forest also returns `variance`, which the inline
# branch dropped and which cannot be recovered post hoc, so patched files have
# no dr_random_forest$variance while newly-run ones do. Nothing in those studies
# reads it - their metrics call only cate_metrics() and hte_test_metrics().
#
# The CI profiles keep tests off; that one IS deliberate (see
# confidence_intervals/README.md).
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
#' @param sl_lib SuperLearner library; NULL skips the SuperLearner arm
#' @param fmla_info oracle formula + parameters; NULL skips the oracle arm
#' @param family gaussian() or binomial(); controls the SuperLearner outcome model
#' @param oracle_link "logit" if fmla_info$fmla is a linear predictor and plogis
#'   must be applied here, "identity" if the formula already includes the link
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
#'   parameter existed - added for continuous/cts_profile.R and friends, not
#'   meant to change production output.
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
                         family = gaussian(), oracle_link = c("identity", "logit"),
                         ipw = NULL, ci = NULL,
                         profile = c("base", "ci", "missing", "ci_mi"),
                         num.threads = NULL, verbose_timing = FALSE,
                         Z_query = NULL) {

  oracle_link <- match.arg(oracle_link)

  profile <- match.arg(profile)
  p <- PROFILES[[profile]]

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
                    oracle_link = oracle_link, tests = p$tests, num.threads = num.threads,
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

# maintainer-endorsed (unsupported) workaround for the fact that grf has no public
# API for an OOB prediction at a counterfactual covariate row: predict(forest)
# without newdata re-reads object$X.orig fresh each call and restricts each row to
# its out-of-bag trees, so pointing X.orig at a perturbed matrix and clearing the
# cached predictions borrows grf's own OOB routine at that perturbed point.
# "Hacky," per grf-labs/grf#307 - may break on a future grf version. R's
# copy-on-modify semantics mean this cannot leak into the caller's forest object;
# only the local copy inside this function is touched. Ported unchanged from
# crossfitting/cf_models.R, where it backs the oob_oob_s arm this file now uses
# as the production default.
oob_predict_counterfactual <- function(forest, X_counterfactual) {
  forest$X.orig <- X_counterfactual
  forest$predictions <- NULL
  forest$debiased.error <- NULL
  predict(forest)$predictions
}

#' Whole-sample OOB nuisance estimation with a regression forest (S-learner)
#'
#' No sample splitting: one S-learner forest on cbind(W, X) supplies
#' Y0.hat/Y1.hat via oob_predict_counterfactual, and two more whole-sample
#' forests supply W.hat and Y.hat.cf, both taken out-of-bag. This is the
#' "oob_oob_s" arm validated in crossfitting/cf_models.R against a T-learner
#' control and a from-scratch manual-API reimplementation
#' (run_all_crossfit_variants there), ported here as the production default
#' in place of the double-crossfit this function used to do.
nuisance_rf <- function(X, Y, W, ipw = NULL, num.threads = NULL) {

  forest <- regression_forest(cbind(W = W, X), Y, sample.weights = ipw,
                              num.threads = num.threads)
  Y0.hat <- oob_predict_counterfactual(forest, cbind(W = 0, X))
  Y1.hat <- oob_predict_counterfactual(forest, cbind(W = 1, X))

  W.hat <- trim_ps(predict(regression_forest(X, W, sample.weights = ipw,
                                             num.threads = num.threads))$predictions)
  Y.hat.cf <- predict(regression_forest(X, Y, sample.weights = ipw,
                                        num.threads = num.threads))$predictions

  Y.hat <- W * Y1.hat + (1 - W) * Y0.hat
  po <- dr_pseudo(Y, W, Y1.hat, Y0.hat, W.hat)

  list(po = po, Y.hat = Y.hat, Y.hat.cf = Y.hat.cf, Y0.hat = Y0.hat, W.hat = W.hat)
}

#' Single leave-one-fold-out nuisance estimation with SuperLearner
#'
#' One split, shared with the stage-2 regression via the same fold_indices
#' (see run_dr_superlearner / stage_2_sl) rather than double-crossfit
#' nuisances feeding a separately-split stage 2 - the "scf_scf" arm validated
#' in crossfitting/cf_models.R.
nuisance_sl <- function(X, Y, W, fold_indices, sl_lib, ipw = NULL,
                        family = gaussian()) {

  binom <- is_binomial(family)

  cross_fits <- future_map(unique(fold_indices), function(fold) {
    in_train <- fold_indices != fold
    in_test <- !in_train

    X_train <- X[in_train, ]
    X_W_train <- cbind(W = W[in_train], X_train)

    Y_lib <- pretest_superlearner(Y[in_train], X_W_train, sl_lib,
                                  if (binom) binomial() else gaussian())
    Y.hat.model <- if (binom) {
      SuperLearner(Y = Y[in_train], X = X_W_train, SL.library = Y_lib,
                   family = binomial(), method = "method.NNloglik",
                   obsWeights = wts(ipw, in_train))
    } else {
      SuperLearner(Y = Y[in_train], X = X_W_train, SL.library = Y_lib,
                   obsWeights = wts(ipw, in_train))
    }

    W_lib <- pretest_superlearner(W[in_train], X_train, sl_lib, binomial())
    W.hat.model <- SuperLearner(W[in_train], X_train, family = binomial(),
                                SL.library = W_lib, method = "method.NNloglik",
                                obsWeights = wts(ipw, in_train))

    X_test <- X[in_test, ]
    Y0.hat <- predict(Y.hat.model, newdata = cbind(W = 0, X_test))$pred
    Y1.hat <- predict(Y.hat.model, newdata = cbind(W = 1, X_test))$pred
    W.hat <- predict(W.hat.model, newdata = X_test)$pred

    # failsafes if SuperLearner returns all-zero predictions
    if (all(Y0.hat == 0) && all(Y1.hat == 0)) {
      warning("SuperLearner failed for Y.hat. Using mean(Y).")
      Y0.hat <- rep(mean(Y[in_train][W[in_train] == 0], na.rm = TRUE), sum(in_test))
      Y1.hat <- rep(mean(Y[in_train][W[in_train] == 1], na.rm = TRUE), sum(in_test))
    }
    if (all(W.hat == 0)) {
      warning("SuperLearner failed for W.hat. Using mean(W).")
      W.hat <- rep(mean(W[in_train], na.rm = TRUE), sum(in_test))
    }

    # clamp propensities away from 0/1, same as nuisance_rf's W.hat
    W.hat <- trim_ps(W.hat)

    W_test <- W[in_test]
    Y.hat <- W_test * Y1.hat + (1 - W_test) * Y0.hat
    po <- dr_pseudo(Y[in_test], W_test, Y1.hat, Y0.hat, W.hat)

    list(fold = fold, po = po, Y.hat = Y.hat, Y0.hat = Y0.hat, W.hat = W.hat)
  }, .options = furrr_options(seed = TRUE))

  scatter_folds(cross_fits, fold_indices, c("po", "Y.hat", "Y0.hat", "W.hat"))
}

#' Drop SuperLearner algorithms that error, warn, or return NA on this data
#'
#' Fits each candidate on its own with a 2-fold inner CV and keeps the survivors.
pretest_superlearner <- function(Y, X, SL.library, family) {
  working_lib <- character()
  removed_lib <- character()
  for (alg in SL.library) {
    fit <- tryCatch(
      SuperLearner(Y = Y, X = X, SL.library = alg, family = family,
                   cvControl = list(V = 2)),
      error = function(e) NULL,
      warning = function(w) NULL
    )
    preds <- if (!is.null(fit) && !is.null(fit$SL.predict)) fit$SL.predict else NULL
    if (!is.null(preds) && any(!is.na(preds))) {
      working_lib <- c(working_lib, alg)
    } else {
      removed_lib <- c(removed_lib, alg)
    }
  }
  if (length(removed_lib) > 0) {
    cat("Removed libraries due to NA/error:\n")
    print(removed_lib)
  }
  if (length(working_lib) == 0) {
    # bug K: every candidate warned/errored on this fold - falling through with
    # character(0) sends an empty SL.library into the caller's live SuperLearner()
    # call, which crashes building a 0-column predictions data.frame() ("arguments
    # imply differing number of rows: 0, 1"). SL.mean is asserted directly, not
    # re-run through the loop above, so it can't recursively trigger the same
    # emptying it's meant to prevent.
    warning("pretest_superlearner: every candidate failed/warned on this fold; ",
            "falling back to SL.mean.")
    working_lib <- "SL.mean"
  }
  working_lib
}

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
stage_2_sl <- function(X, po, fold_indices, fold_list, sl_lib, ipw = NULL) {
  n_obs <- nrow(X)
  single <- is.vector(po)

  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold <- !in_train

    y_train <- if (single) po[in_train] else po[in_train, fold]
    po_lib <- pretest_superlearner(y_train, X[in_train, ], sl_lib, gaussian())

    lib <- po_lib

    po_model <- SuperLearner(y_train, X[in_train, ], family = gaussian(),
                             SL.library = lib, obsWeights = wts(ipw, in_train))
    list(fold = fold, predictions = predict(po_model, newdata = X[in_fold, ])$pred)
  }, .options = furrr_options(seed = TRUE))

  tau <- rep(NA, n_obs)
  for (result in tau_results) tau[fold_indices == result$fold] <- result$predictions
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
run_dr_oracle <- function(X, Y, W, fmla_info, ipw = NULL,
                          oracle_link = c("identity", "logit"), tests = TRUE,
                          num.threads = NULL, Z_query = NULL) {
  n_obs <- nrow(X)
  # "logit" means the formula is a linear predictor and plogis belongs here;
  # "identity" means the formula already returns the outcome mean
  link <- if (match.arg(oracle_link) == "logit") plogis else identity

  X <- as.data.frame(X)
  list2env(fmla_info$params, envir = environment())
  fmla <- parse(text = fmla_info$fmla)

  W_temp <- rep(1, n_obs)
  Y1.hat <- link(eval(fmla, envir = list2env(c(list(W = W_temp), X))))

  W_temp <- rep(0, n_obs)
  Y0.hat <- link(eval(fmla, envir = list2env(c(list(W = W_temp), X))))

  Y.hat <- link(eval(fmla, envir = list2env(c(list(W = W), X))))
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

  forest <- regression_forest(cbind(W = W, X), Y, sample.weights = ipw,
                              num.threads = num.threads)
  Y0.hat <- oob_predict_counterfactual(forest, cbind(W = 0, X))
  Y1.hat <- oob_predict_counterfactual(forest, cbind(W = 1, X))

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

  out <- list(tau = tau)
  if (tests) {
    out$BLP_whole <- run_blp_whole(Y, W, nuisances$W.hat, nuisances$Y0.hat, tau)
    out$independence_cate <- run_independence_test_whole(X, tau)
    out$independence_po <- run_independence_test_whole(X, nuisances$po)
  }
  out
}

# ---- post-estimation heterogeneity tests ------------------------------------

# Best Linear Predictor of the CATE (GenericML). Row 4, column 2 of the returned
# coefficient block is the p-value the metrics scripts read.
#
# bug L: GenericML::BLP() regresses on beta.2 = (W - W.hat) * (tau - mean(tau)).
# When tau is exactly constant (a degenerate/near-constant CATE fit - seen with
# scenario 9's low-amplitude cos(X4) effect at small n, especially once
# pretest_superlearner has whittled a fold's library down to 1-2 survivors),
# beta.2 becomes identically zero, lm() marks it aliased (coef = NA), and
# sandwich::vcovHC() drops that coefficient's row/column entirely rather than
# keeping it as NA - so GenericML's internal indexing by name throws "subscript
# out of bounds". NULL is a deliberate fallback, not just a safe default:
# hte_test_metrics() (R/metrics.R) already maps BLP_whole = NULL to BLP_p = NA,
# and NA is the statistically correct answer here - beta.2 has no fitted
# coefficient to attach a p-value to when tau has zero variance.
run_blp_whole <- function(Y, W, W.hat, Y0.hat, tau) {
  tryCatch(
    BLP(Y, W, W.hat, Y0.hat, tau)$coefficients[, c(1, 4)],
    error = function(e) {
      warning("run_blp_whole: BLP() failed (likely a constant/degenerate tau); ",
              "returning NULL. ", conditionMessage(e))
      NULL
    }
  )
}

# Omnibus independence test of the estimated CATEs against the covariates
run_independence_test_whole <- function(X, tau) {
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
      method = "independence_test"
    )
  }, error = function(e) {
    list(p_value = 1, statistic = 0, method = "independence_test_failed")
  })
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
