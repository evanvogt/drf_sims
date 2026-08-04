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

# TEMPORARY (bug F). stage_2_sl computes a pretested library and then, in the
# matrix branch, passes the untested one. Every caller passes a matrix, so the
# pretested library has never been used in stage 2 in any study, and the vector
# branch below is dead code. FALSE reproduces that; Step 8 flips it to TRUE and
# deletes the option.
PRETEST_STAGE2 <- TRUE

# ---- orchestration profiles -------------------------------------------------
#
#                        base      ci        missing
#  causal forest variance  no       yes       yes
#  causal forest tests     yes      no        yes
#  dr_random_forest tests  yes      no        no      <- divergence, see NOTE
#  oracle / semi tests     yes      no        yes
#  SuperLearner arm        yes      no        yes (skipped if X still has NAs)
#  half-sample bootstrap   no       yes       no
#  nuisance row means      yes      no        yes
#
# NOTE: in the missing-data and CI variants dr_random_forest is built inline as
# list(tau = stage_2_rf(...)) and so carries no BLP or independence test, while
# the base studies call run_dr_random_forest and get both. The metrics scripts
# then record BLP_p as NA for that one model in those studies. Nothing marks this
# as deliberate and it looks like copy-paste drift, but it is preserved here
# because changing it would move published numbers. Worth a decision.

PROFILES <- list(
  base    = list(cf_variance = FALSE, tests = TRUE,  dr_rf_tests = TRUE,
                 aggregate_nuisances = TRUE),
  ci      = list(cf_variance = TRUE,  tests = FALSE, dr_rf_tests = FALSE,
                 aggregate_nuisances = FALSE),
  missing = list(cf_variance = TRUE,  tests = TRUE,  dr_rf_tests = FALSE,
                 aggregate_nuisances = TRUE),
  # missing/ci_example: as "ci", but it kept the nuisance row means
  ci_mi   = list(cf_variance = TRUE,  tests = FALSE, dr_rf_tests = FALSE,
                 aggregate_nuisances = TRUE)
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
#' @param n_folds number of crossfitting folds V
#' @param sl_lib SuperLearner library; NULL skips the SuperLearner arm
#' @param fmla_info oracle formula + parameters; NULL skips the oracle arm
#' @param family gaussian() or binomial(); controls the SuperLearner outcome model
#' @param oracle_link "logit" if fmla_info$fmla is a linear predictor and plogis
#'   must be applied here, "identity" if the formula already includes the link
#' @param ipw optional length-n weights for the missing-data IPW arm
#' @param ci NULL, or list(boot = , sf = , alpha = ) to add half-sample bootstrap CIs
#' @param profile "base", "ci" or "missing" - see PROFILES
#'
#' Each study's *_models.R is a thin shim defining `run_all_cate_methods` with
#' that study's historical signature and forwarding to this. The shared function
#' is named differently so those shims do not recurse into themselves.
cate_methods <- function(data, n_folds = 10, sl_lib = NULL, fmla_info = NULL,
                         family = gaussian(), oracle_link = c("identity", "logit"),
                         ipw = NULL, ci = NULL,
                         profile = c("base", "ci", "missing", "ci_mi")) {

  oracle_link <- match.arg(oracle_link)

  profile <- match.arg(profile)
  p <- PROFILES[[profile]]

  X <- as.matrix(data[, -c(1:2)])
  Y <- data$Y
  W <- data$W
  n_obs <- nrow(X)

  fold_indices <- sort(seq(n_obs) %% n_folds) + 1
  fold_list <- unique(fold_indices)
  fold_pairs <- utils::combn(fold_list, 2, simplify = FALSE)

  results <- list()

  cat("Computing nuisance functions...\n")
  nuisances_rf <- nuisance_rf(X, Y, W, fold_indices, fold_pairs, ipw,
                              aggregate = p$aggregate_nuisances)

  cat("Running Causal Forest...\n")
  results$causal_forest <- run_causal_forest(X, Y, W, nuisances_rf, fold_indices,
                                             fold_list, ipw,
                                             variance = p$cf_variance,
                                             tests = p$tests)
  if (!is.null(ci)) {
    cat("Running Causal Forest bootstrap... \n")
    results$causal_forest <- c(results$causal_forest,
      cf_half_boot(X, Y, W, nuisances_rf, results$causal_forest$tau,
                   ci$boot, ci$sf, ci$alpha, fold_indices, fold_list))
  }

  cat("Running DR Random Forest...\n")
  results$dr_random_forest <- if (p$dr_rf_tests) {
    run_dr_random_forest(X, Y, W, nuisances_rf, fold_indices, fold_list, ipw)
  } else {
    list(tau = stage_2_rf(X, nuisances_rf$po_matrix, fold_indices, fold_list, ipw))
  }
  if (!is.null(ci)) {
    cat("Running DR RF bootstrap... \n")
    results$dr_random_forest <- c(results$dr_random_forest,
      rf_half_boot(X, Y, W, nuisances_rf$po_matrix, results$dr_random_forest$tau,
                   ci$boot, ci$sf, ci$alpha, fold_indices, fold_list))
  }

  # the missing-data variant skips the arms that need a complete covariate matrix
  complete_X <- !anyNA(X)

  if (!is.null(fmla_info) && (profile != "missing" || complete_X)) {
    cat("Running DR Oracle...\n")
    results$dr_oracle <- run_dr_oracle(X, Y, W, fmla_info, fold_indices, fold_list,
                                       ipw, oracle_link = oracle_link,
                                       tests = p$tests)
    if (!is.null(ci)) {
      cat("Runnings Oracle bootstrap...\n")
      results$dr_oracle <- c(results$dr_oracle,
        rf_half_boot(X, Y, W, results$dr_oracle$po, results$dr_oracle$tau,
                     ci$boot, ci$sf, ci$alpha, fold_indices, fold_list))
    }
  }

  if (profile != "missing" || complete_X) {
    cat("Running DR Semi-Oracle...\n")
    results$dr_semi_oracle <- run_dr_semi_oracle(X, Y, W, fold_indices, fold_list,
                                                 ipw, tests = p$tests)
    if (!is.null(ci)) {
      cat("Running Semi-Oracle bootstrap...\n")
      results$dr_semi_oracle <- c(results$dr_semi_oracle,
        rf_half_boot(X, Y, W, results$dr_semi_oracle$po, results$dr_semi_oracle$tau,
                     ci$boot, ci$sf, ci$alpha, fold_indices, fold_list))
    }
  }

  if (!is.null(sl_lib) && (profile != "missing" || complete_X)) {
    cat("Running DR SuperLearner...\n")
    X <- as.data.frame(X)
    nuisances_sl <- nuisance_sl(X, Y, W, fold_indices, fold_pairs, sl_lib, ipw,
                                family = family)
    results$dr_superlearner <- run_dr_superlearner(X, Y, W, nuisances_sl,
                                                   fold_indices, fold_list,
                                                   fold_pairs, sl_lib, ipw,
                                                   tests = p$tests)
    results$nuisances_sl <- nuisances_sl
  }

  results$nuisances_rf <- nuisances_rf
  results$fold_indices <- fold_indices

  results
}

# ---- stage 1: nuisance estimation -------------------------------------------

#' Double-crossfit nuisance estimation with regression forests
#'
#' Fits over all C(V, 2) fold pairs, so column k of each returned matrix is
#' untouched by fold k. `aggregate` adds the row-mean summaries the BLP and
#' independence tests consume; the CI variant never ran those and so never
#' carried them in its saved output.
nuisance_rf <- function(X, Y, W, fold_indices, fold_pairs, ipw = NULL,
                        aggregate = TRUE) {

  cross_fits <- future_map(seq_along(fold_pairs), function(i) {
    fold_pair <- fold_pairs[[i]]
    in_train <- !(fold_indices %in% fold_pair)
    in_test <- !in_train

    Y.hat.model <- regression_forest(cbind(W[in_train], X[in_train, ]), Y[in_train],
                                     sample.weights = wts(ipw, in_train))
    Y.hat.cf.model <- regression_forest(X[in_train, ], Y[in_train],
                                        sample.weights = wts(ipw, in_train))
    W.hat.model <- regression_forest(X[in_train, ], W[in_train],
                                     sample.weights = wts(ipw, in_train))

    X_test <- X[in_test, ]

    Y0.hat <- predict(Y.hat.model, newdata = cbind(W = 0, X_test))$predictions
    Y1.hat <- predict(Y.hat.model, newdata = cbind(W = 1, X_test))$predictions
    Y.hat.cf <- predict(Y.hat.cf.model, newdata = X_test)$predictions
    W.hat <- predict(W.hat.model, newdata = X_test)$predictions

    W_test <- W[in_test]
    Y.hat <- W_test * Y1.hat + (1 - W_test) * Y0.hat
    po <- dr_pseudo(Y[in_test], W_test, Y1.hat, Y0.hat, W.hat)

    list(po = po, Y.hat = Y.hat, Y0.hat = Y0.hat, Y.hat.cf = Y.hat.cf,
         W.hat = W.hat, fold_pair = fold_pair)
  }, .options = furrr_options(seed = TRUE))

  fold_list <- unique(fold_indices)
  mats <- lapply(c("po", "Y.hat", "Y.hat.cf", "Y0.hat", "W.hat"), function(nm) {
    collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, nm)
  })
  names(mats) <- c("po_matrix", "Y.hat_matrix", "Y.hat.cf_matrix",
                   "Y0.hat_matrix", "W.hat_matrix")

  if (!aggregate) {
    # field order as the CI variant emitted it
    return(mats[c("po_matrix", "Y.hat_matrix", "Y.hat.cf_matrix",
                  "Y0.hat_matrix", "W.hat_matrix")])
  }

  list(
    po_matrix = mats$po_matrix,
    po = rowMeans(mats$po_matrix, na.rm = TRUE),
    Y.hat_matrix = mats$Y.hat_matrix,
    Y.hat = rowMeans(mats$Y.hat_matrix, na.rm = TRUE),
    Y.hat.cf_matrix = mats$Y.hat.cf_matrix,
    Y.hat.cf = rowMeans(mats$Y.hat.cf_matrix, na.rm = TRUE),
    Y0.hat_matrix = mats$Y0.hat_matrix,
    Y0.hat = rowMeans(mats$Y0.hat_matrix, na.rm = TRUE),
    W.hat_matrix = mats$W.hat_matrix,
    W.hat = rowMeans(mats$W.hat_matrix, na.rm = TRUE)
  )
}

#' Double-crossfit nuisance estimation with SuperLearner
nuisance_sl <- function(X, Y, W, fold_indices, fold_pairs, sl_lib, ipw = NULL,
                        family = gaussian()) {

  binom <- is_binomial(family)

  cross_fits <- future_map(seq_along(fold_pairs), function(i) {
    fold_pair <- fold_pairs[[i]]
    in_train <- !(fold_indices %in% fold_pair)
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

    # trim extreme propensities. NOTE: done only on this SuperLearner path, not
    # in nuisance_rf - an asymmetry inherited from the original code. It is a
    # no-op while W is randomised 0.5, but it is not obviously intended.
    W.hat[W.hat < 0.05] <- 0.05
    W.hat[W.hat > 0.95] <- 0.95

    W_test <- W[in_test]
    Y.hat <- W_test * Y1.hat + (1 - W_test) * Y0.hat
    po <- dr_pseudo(Y[in_test], W_test, Y1.hat, Y0.hat, W.hat)

    list(po = po, Y.hat = Y.hat, Y0.hat = Y0.hat, W.hat = W.hat,
         fold_pair = fold_pair)
  }, .options = furrr_options(seed = TRUE))

  fold_list <- unique(fold_indices)
  po_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "po")
  Y.hat_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "Y.hat")
  Y0.hat_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "Y0.hat")
  W.hat_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "W.hat")

  list(
    po_matrix = po_matrix,
    po = rowMeans(po_matrix, na.rm = TRUE),
    Y.hat_matrix = Y.hat_matrix,
    Y.hat = rowMeans(Y.hat_matrix, na.rm = TRUE),
    Y0.hat_matrix = Y0.hat_matrix,
    Y0.hat = rowMeans(Y0.hat_matrix, na.rm = TRUE),
    W.hat_matrix = W.hat_matrix,
    W.hat = rowMeans(W.hat_matrix, na.rm = TRUE)
  )
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
  working_lib
}

# ---- stage 2: final CATE regression -----------------------------------------

#' Crossfit second stage with a regression forest
#'
#' `po` is either the n x V double-crossfitting matrix (column k valid for fold k)
#' or a plain n-vector, as produced by the oracle arms.
stage_2_rf <- function(X, po, fold_indices, fold_list, ipw = NULL) {
  n_obs <- nrow(X)
  single <- is.vector(po)

  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold <- !in_train

    y_train <- if (single) po[in_train] else po[in_train, fold]
    forest <- regression_forest(X[in_train, ], y_train,
                                sample.weights = wts(ipw, in_train))
    list(fold = fold, predictions = predict(forest, newdata = X[in_fold, ])$predictions)
  }, .options = furrr_options(seed = TRUE))

  tau <- rep(NA, n_obs)
  for (result in tau_results) tau[fold_indices == result$fold] <- result$predictions
  tau
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

    # bug F: the matrix branch historically passed the untested sl_lib even though
    # po_lib had just been computed. Every caller passes a matrix, so the pretested
    # library has never actually been used here. PRETEST_STAGE2 gates the fix.
    lib <- if (single || PRETEST_STAGE2) po_lib else sl_lib

    po_model <- SuperLearner(y_train, X[in_train, ], family = gaussian(),
                             SL.library = lib, obsWeights = wts(ipw, in_train))
    list(fold = fold, predictions = predict(po_model, newdata = X[in_fold, ])$pred)
  }, .options = furrr_options(seed = TRUE))

  tau <- rep(NA, n_obs)
  for (result in tau_results) tau[fold_indices == result$fold] <- result$predictions
  tau
}

# ---- the estimators ---------------------------------------------------------

#' Fold-wise causal forest using pre-computed double-crossfit nuisances
run_causal_forest <- function(X, Y, W, nuisances, fold_indices, fold_list,
                              ipw = NULL, variance = FALSE, tests = TRUE) {
  n_obs <- nrow(X)

  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold <- !in_train

    forest <- causal_forest(X[in_train, ], Y[in_train], W[in_train],
                            nuisances$Y.hat.cf_matrix[in_train, fold],
                            nuisances$W.hat_matrix[in_train, fold],
                            sample.weights = wts(ipw, in_train))

    pred <- predict(forest, newdata = X[in_fold, ], estimate.variance = variance)
    list(fold = fold, tau = pred$predictions,
         variance = if (variance) pred$variance.estimates else NULL)
  }, .options = furrr_options(seed = TRUE))

  tau <- rep(NA, n_obs)
  tau_var <- rep(NA, n_obs)
  for (result in tau_results) {
    in_fold <- fold_indices == result$fold
    tau[in_fold] <- result$tau
    if (variance) tau_var[in_fold] <- result$variance
  }

  out <- list(tau = tau)
  if (variance) out$variance <- tau_var
  if (tests) {
    out$BLP_whole <- run_blp_whole(Y, W, nuisances$W.hat, nuisances$Y0.hat, tau)
    out$independence_cate <- run_independence_test_whole(X, tau)
    out$independence_po <- run_independence_test_whole(X, nuisances$po)
  }
  out
}

#' DR-learner with a regression-forest second stage
run_dr_random_forest <- function(X, Y, W, nuisances, fold_indices, fold_list,
                                 ipw = NULL, tests = TRUE) {
  tau <- stage_2_rf(X, nuisances$po_matrix, fold_indices, fold_list, ipw)

  out <- list(tau = tau)
  if (tests) {
    out$BLP_whole <- run_blp_whole(Y, W, nuisances$W.hat, nuisances$Y0.hat, tau)
    out$independence_cate <- run_independence_test_whole(X, tau)
    out$independence_po <- run_independence_test_whole(X, nuisances$po)
  }
  out
}

#' DR-learner with the true outcome model and a known propensity of 0.5
run_dr_oracle <- function(X, Y, W, fmla_info, fold_indices, fold_list, ipw = NULL,
                          oracle_link = c("identity", "logit"), tests = TRUE) {
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
  tau <- stage_2_rf(X, po, fold_indices, fold_list, ipw)

  out <- list(tau = tau, po = po, Y0.hat = Y0.hat)
  if (tests) {
    out$BLP_whole <- run_blp_whole(Y, W, W.hat, Y0.hat, tau)
    out$independence_cate <- run_independence_test_whole(X, tau)
    out$independence_po <- run_independence_test_whole(X, po)
  }
  out
}

#' DR-learner with a known propensity of 0.5 but an estimated outcome model
run_dr_semi_oracle <- function(X, Y, W, fold_indices, fold_list, ipw = NULL,
                               tests = TRUE) {
  n_obs <- nrow(X)
  W.hat <- rep(0.5, n_obs)

  cross_fits <- future_map(seq_len(length(fold_list)), function(fold) {
    in_train <- !(fold_indices == fold)
    in_test <- which(!in_train)

    Y.hat.model <- regression_forest(cbind(W[in_train], X[in_train, ]), Y[in_train],
                                     sample.weights = wts(ipw, in_train))

    X_test <- X[in_test, ]
    Y0.hat <- predict(Y.hat.model, newdata = cbind(W = 0, X_test))$predictions
    Y1.hat <- predict(Y.hat.model, newdata = cbind(W = 1, X_test))$predictions

    W_test <- W[in_test]
    po <- dr_pseudo(Y[in_test], W_test, Y1.hat, Y0.hat, W.hat)

    list(po = po, Y0.hat = Y0.hat, fold = fold)
  }, .options = furrr_options(seed = TRUE))

  po <- rep(NA, n_obs)
  Y0.hat <- rep(NA, n_obs)
  for (i in seq_along(cross_fits)) {
    in_fold <- fold_indices == cross_fits[[i]]$fold
    po[in_fold] <- cross_fits[[i]]$po
    Y0.hat[in_fold] <- cross_fits[[i]]$Y0.hat
  }

  tau <- stage_2_rf(X, po, fold_indices, fold_list, ipw)

  out <- list(tau = tau, po = po, Y0.hat = Y0.hat)
  if (tests) {
    out$BLP_whole <- run_blp_whole(Y, W, W.hat, Y0.hat, tau)
    out$independence_cate <- run_independence_test_whole(X, tau)
    out$independence_po <- run_independence_test_whole(X, po)
  }
  out
}

#' DR-learner with SuperLearner nuisances and second stage
run_dr_superlearner <- function(X, Y, W, nuisances, fold_indices, fold_list,
                                fold_pairs, sl_lib, ipw = NULL, tests = TRUE) {
  tau <- stage_2_sl(X, nuisances$po_matrix, fold_indices, fold_list, sl_lib, ipw)

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
run_blp_whole <- function(Y, W, W.hat, Y0.hat, tau) {
  BLP(Y, W, W.hat, Y0.hat, tau)$coefficients[, c(1, 4)]
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
