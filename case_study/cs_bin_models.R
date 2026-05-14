##########
# Title: multi-arm CATE models for binary outcome - factorial platform trial case study
# W must be a factor variable with K >= 2 levels; first level is the reference arm
# Returns results for ALL pairwise comparisons
##########
library(coin)
library(dplyr)
library(furrr)
library(GenericML)
library(grf)
library(SuperLearner)

# ---- Utilities (self-contained) ----

collate_predictions <- function(fold_list, fold_pairs, fold_indices, reslist, target) {
  lapply(fold_list, function(fold) {
    predictions <- rep(NA, length(fold_indices))
    for (j in seq_along(fold_pairs)) {
      if (fold %in% fold_pairs[[j]]) {
        predictions[fold_indices %in% fold_pairs[[j]]] <- reslist[[j]][[target]]
      }
    }
    predictions[fold_indices == fold] <- NA
    predictions
  }) %>% simplify2array()
}

stage_2_rf <- function(X, po, fold_indices, fold_list) {
  n_obs  <- nrow(X)
  single <- is.vector(po)

  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold     <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold  <- !in_train
    forest   <- if (single) regression_forest(X[in_train, ], po[in_train])
                else        regression_forest(X[in_train, ], po[in_train, fold])
    list(fold = fold, predictions = predict(forest, newdata = X[in_fold, ])$predictions)
  }, .options = furrr_options(seed = TRUE))

  tau <- rep(NA_real_, n_obs)
  for (r in tau_results) tau[fold_indices == r$fold] <- r$predictions
  tau
}

stage_2_sl <- function(X, po, fold_indices, fold_list, sl_lib) {
  n_obs  <- nrow(X)
  single <- is.vector(po)

  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold     <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold  <- !in_train
    po_train <- if (single) po[in_train] else po[in_train, fold]
    po_lib   <- pretest_superlearner(po_train, X[in_train, ], sl_lib, gaussian())
    model    <- SuperLearner(po_train, X[in_train, ], family = gaussian(), SL.library = po_lib)
    list(fold = fold, predictions = predict(model, newdata = X[in_fold, ])$pred)
  }, .options = furrr_options(seed = TRUE))

  tau <- rep(NA_real_, n_obs)
  for (r in tau_results) tau[fold_indices == r$fold] <- r$predictions
  tau
}

pretest_superlearner <- function(Y, X, SL.library, family) {
  working_lib <- character()
  removed_lib <- character()
  for (alg in SL.library) {
    fit <- tryCatch(
      SuperLearner(Y = Y, X = X, SL.library = alg, family = family, cvControl = list(V = 2)),
      error = function(e) NULL, warning = function(w) NULL
    )
    preds <- if (!is.null(fit) && !is.null(fit$SL.predict)) fit$SL.predict else NULL
    if (!is.null(preds) && any(!is.na(preds))) working_lib <- c(working_lib, alg)
    else                                        removed_lib <- c(removed_lib, alg)
  }
  if (length(removed_lib) > 0) { cat("Removed SL libraries:\n"); print(removed_lib) }
  working_lib
}

run_blp_whole <- function(Y, W, W.hat, Y0.hat, tau) {
  BLP(Y, W, W.hat, Y0.hat, tau)$coefficients[, c(1, 4)]
}

run_independence_test_whole <- function(X, tau) {
  tryCatch({
    res <- coin::independence_test(tau ~ ., data = data.frame(tau = tau, X), teststat = "quadratic")
    list(p_value = coin::pvalue(res), statistic = coin::statistic(res), method = "independence_test")
  }, error = function(e) {
    list(p_value = 1, statistic = 0, method = "independence_test_failed")
  })
}

# ---- Helpers ----

# Sign convention: pair_name(arm_k, arm_j) = tau(arm_k) - tau(arm_j)
pair_name <- function(arm_k, arm_j) paste0(arm_k, "_vs_", arm_j)

# Derive all-pairwise contrasts from multi_arm_causal_forest vs-reference predictions
# macf_preds: n x (K-1) matrix, cols named e.g. "A2 - A1", "A3 - A1"
derive_all_pairwise <- function(macf_preds, arm_levels, all_pairs) {
  if (!is.matrix(macf_preds)) {
    macf_preds <- matrix(macf_preds, ncol = 1,
                         dimnames = list(NULL, paste0(arm_levels[2], " - ", arm_levels[1])))
  }

  ref_arm     <- arm_levels[1]
  active_arms <- arm_levels[-1]

  vs_ref <- setNames(
    lapply(active_arms, function(arm) macf_preds[, paste0(arm, " - ", ref_arm)]),
    active_arms
  )

  setNames(
    lapply(all_pairs, function(pair) {
      arm_k <- pair[1]; arm_j <- pair[2]
      if      (arm_k == ref_arm) -vs_ref[[arm_j]]
      else if (arm_j == ref_arm)  vs_ref[[arm_k]]
      else                        vs_ref[[arm_k]] - vs_ref[[arm_j]]
    }),
    sapply(all_pairs, function(p) pair_name(p[1], p[2]))
  )
}

# ---- Main orchestrator ----

run_all_cate_methods_multiarm <- function(data, outcome = "Y", n_folds = 10, sl_lib = NULL) {
  Y <- data[[outcome]]
  W <- data$W
  X <- as.matrix(data[, !names(data) %in% c(outcome, "W")])
  n_obs <- nrow(X)

  if (!is.factor(W)) stop("W must be a factor variable with named levels")

  arm_levels <- levels(W)
  n_arms     <- length(arm_levels)
  all_pairs  <- combn(arm_levels, 2, simplify = FALSE)

  fold_indices <- sample(sort(seq(n_obs) %% n_folds) + 1)
  fold_list    <- 1:n_folds
  fold_pairs   <- utils::combn(fold_list, 2, simplify = FALSE)

  results <- list()

  message("Computing multi-arm nuisance functions (RF)...")
  nuisances_rf <- nuisance_rf_multiarm(X, Y, W, fold_indices, fold_pairs, arm_levels, all_pairs)

  message("Running Multi-arm Causal Forest...")
  results$causal_forest <- run_causal_forest_multiarm(X, Y, W, nuisances_rf, fold_indices, fold_list,
                                                       arm_levels, all_pairs)

  message("Running Multi-arm DR Random Forest...")
  results$dr_random_forest <- run_dr_random_forest_multiarm(X, Y, W, nuisances_rf, fold_indices,
                                                             fold_list, arm_levels, all_pairs)

  if (!is.null(sl_lib)) {
    message("Running Multi-arm DR SuperLearner...")
    X_df         <- as.data.frame(X)
    nuisances_sl <- nuisance_sl_multiarm(X_df, Y, W, fold_indices, fold_pairs, arm_levels, all_pairs, sl_lib)
    results$dr_superlearner <- run_dr_superlearner_multiarm(X_df, Y, W, nuisances_sl, fold_indices,
                                                             fold_list, arm_levels, all_pairs, sl_lib)
    results$nuisances_sl <- nuisances_sl
  }

  results$nuisances_rf <- nuisances_rf
  results$fold_indices <- fold_indices
  results$arm_levels   <- arm_levels
  results$all_pairs    <- all_pairs

  return(results)
}

# ---- Method functions ----

# Double-crossfit RF nuisances for multi-arm binary outcome
nuisance_rf_multiarm <- function(X, Y, W, fold_indices, fold_pairs, arm_levels, all_pairs) {
  n_obs      <- nrow(X)
  n_arms     <- length(arm_levels)
  pair_names <- sapply(all_pairs, function(p) pair_name(p[1], p[2]))
  fold_list  <- sort(unique(fold_indices))

  W_dummies <- model.matrix(~ W - 1)
  colnames(W_dummies) <- arm_levels

  cross_fits <- future_map(seq_along(fold_pairs), function(i) {
    fold_pair <- fold_pairs[[i]]
    in_train  <- !(fold_indices %in% fold_pair)
    in_test   <- !in_train
    n_test    <- sum(in_test)

    X_train     <- X[in_train, ]
    X_test      <- X[in_test, ]
    W_dum_train <- W_dummies[in_train, , drop = FALSE]
    W_test      <- W[in_test]

    mu_model       <- regression_forest(cbind(W_dum_train, X_train), Y[in_train])
    Y.hat.cf.model <- regression_forest(X_train, Y[in_train])
    Y.hat.cf       <- predict(Y.hat.cf.model, newdata = X_test)$predictions

    pf         <- probability_forest(X_train, W[in_train])
    e_hat_test <- predict(pf, newdata = X_test)$predictions
    if (is.null(colnames(e_hat_test)) || !all(colnames(e_hat_test) == arm_levels)) {
      colnames(e_hat_test) <- arm_levels
    }

    mu_hat_per_arm <- setNames(
      lapply(arm_levels, function(arm) {
        W_pred <- matrix(0, nrow = n_test, ncol = n_arms)
        colnames(W_pred) <- arm_levels
        W_pred[, arm] <- 1
        predict(mu_model, newdata = cbind(W_pred, X_test))$predictions
      }),
      arm_levels
    )

    W_test_char <- as.character(W_test)
    mu_observed <- mapply(function(arm, j) mu_hat_per_arm[[arm]][j],
                          W_test_char, seq_along(W_test_char))

    pos <- setNames(
      lapply(all_pairs, function(pair) {
        arm_k <- pair[1]; arm_j <- pair[2]
        e_k   <- e_hat_test[, arm_k]
        e_j   <- e_hat_test[, arm_j]
        mu_k  <- mu_hat_per_arm[[arm_k]]
        mu_j  <- mu_hat_per_arm[[arm_j]]
        ind_k <- as.integer(W_test == arm_k)
        ind_j <- as.integer(W_test == arm_j)
        cate  <- mu_k - mu_j
        cate + (ind_k / e_k - ind_j / e_j) * (Y[in_test] - mu_observed)
      }),
      pair_names
    )

    e_per_arm <- setNames(lapply(arm_levels, function(arm) e_hat_test[, arm]), arm_levels)

    list(pos = pos, Y.hat.cf = Y.hat.cf, e_per_arm = e_per_arm,
         mu_hat_per_arm = mu_hat_per_arm, fold_pair = fold_pair)
  }, .options = furrr_options(seed = TRUE))

  po_matrices <- setNames(
    lapply(pair_names, function(pname) {
      collate_predictions(fold_list, fold_pairs, fold_indices,
                          lapply(cross_fits, function(cf) list(po = cf$pos[[pname]], fold_pair = cf$fold_pair)),
                          "po")
    }),
    pair_names
  )

  Y.hat.cf_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "Y.hat.cf")

  W.hat <- do.call(cbind, setNames(
    lapply(arm_levels, function(arm) {
      e_mat <- collate_predictions(fold_list, fold_pairs, fold_indices,
                                   lapply(cross_fits, function(cf) list(e_arm = cf$e_per_arm[[arm]], fold_pair = cf$fold_pair)),
                                   "e_arm")
      rowMeans(e_mat, na.rm = TRUE)
    }),
    arm_levels
  ))

  mu_per_arm <- setNames(
    lapply(arm_levels, function(arm) {
      mu_mat <- collate_predictions(fold_list, fold_pairs, fold_indices,
                                    lapply(cross_fits, function(cf) list(mu_arm = cf$mu_hat_per_arm[[arm]], fold_pair = cf$fold_pair)),
                                    "mu_arm")
      rowMeans(mu_mat, na.rm = TRUE)
    }),
    arm_levels
  )

  list(
    po_matrices     = po_matrices,
    po              = setNames(lapply(pair_names, function(pn) rowMeans(po_matrices[[pn]], na.rm = TRUE)), pair_names),
    Y.hat.cf_matrix = Y.hat.cf_matrix,
    Y.hat.cf        = rowMeans(Y.hat.cf_matrix, na.rm = TRUE),
    mu_per_arm      = mu_per_arm,
    W.hat           = W.hat
  )
}

# Multi-arm causal forest via multi_arm_causal_forest + derive_all_pairwise
run_causal_forest_multiarm <- function(X, Y, W, nuisances, fold_indices, fold_list,
                                        arm_levels, all_pairs) {
  n_obs      <- nrow(X)
  pair_names <- sapply(all_pairs, function(p) pair_name(p[1], p[2]))

  fold_results <- future_map(seq_along(fold_list), function(i) {
    fold     <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold  <- !in_train

    forest <- multi_arm_causal_forest(
      X[in_train, ], Y[in_train], W[in_train],
      Y.hat = nuisances$Y.hat.cf_matrix[in_train, i],
      W.hat = nuisances$W.hat[in_train, , drop = FALSE]
    )

    pred     <- predict(forest, newdata = X[in_fold, ], drop = TRUE)$predictions
    pairwise <- derive_all_pairwise(pred, arm_levels, all_pairs)
    list(fold = fold, pairwise = pairwise)
  }, .options = furrr_options(seed = TRUE))

  tau_list <- setNames(lapply(pair_names, function(pn) rep(NA_real_, n_obs)), pair_names)
  for (r in fold_results) {
    in_fold <- fold_indices == r$fold
    for (pn in pair_names) tau_list[[pn]][in_fold] <- r$pairwise[[pn]]
  }

  list(
    tau         = tau_list,
    diagnostics = run_diagnostics_pairwise(X, Y, W, arm_levels, all_pairs, tau_list, nuisances)
  )
}

# DR random forest (multi-arm): stage_2_rf per pair
run_dr_random_forest_multiarm <- function(X, Y, W, nuisances, fold_indices, fold_list,
                                           arm_levels, all_pairs) {
  pair_names <- sapply(all_pairs, function(p) pair_name(p[1], p[2]))

  tau_list <- setNames(
    lapply(pair_names, function(pn) stage_2_rf(X, nuisances$po_matrices[[pn]], fold_indices, fold_list)),
    pair_names
  )

  list(
    tau         = tau_list,
    diagnostics = run_diagnostics_pairwise(X, Y, W, arm_levels, all_pairs, tau_list, nuisances)
  )
}

# Double-crossfit SL nuisances for multi-arm binary outcome
nuisance_sl_multiarm <- function(X, Y, W, fold_indices, fold_pairs, arm_levels, all_pairs, sl_lib) {
  n_obs      <- nrow(X)
  n_arms     <- length(arm_levels)
  pair_names <- sapply(all_pairs, function(p) pair_name(p[1], p[2]))
  fold_list  <- sort(unique(fold_indices))

  W_dummies <- model.matrix(~ W - 1)
  colnames(W_dummies) <- arm_levels

  cross_fits <- future_map(seq_along(fold_pairs), function(i) {
    fold_pair <- fold_pairs[[i]]
    in_train  <- !(fold_indices %in% fold_pair)
    in_test   <- !in_train
    n_test    <- sum(in_test)

    X_train     <- X[in_train, ]
    X_test      <- X[in_test, ]
    W_dum_train <- W_dummies[in_train, , drop = FALSE]
    W_test      <- W[in_test]
    X_W_train   <- cbind(W_dum_train, X_train)

    mu_lib   <- pretest_superlearner(Y[in_train], X_W_train, sl_lib, binomial())
    mu_model <- SuperLearner(Y = Y[in_train], X = X_W_train, SL.library = mu_lib,
                             family = binomial(), method = "method.NNloglik")

    pf         <- probability_forest(X_train, W[in_train])
    e_hat_test <- predict(pf, newdata = X_test)$predictions
    if (is.null(colnames(e_hat_test)) || !all(colnames(e_hat_test) == arm_levels)) {
      colnames(e_hat_test) <- arm_levels
    }

    mu_hat_per_arm <- setNames(
      lapply(arm_levels, function(arm) {
        W_pred <- matrix(0, nrow = n_test, ncol = n_arms)
        colnames(W_pred) <- arm_levels
        W_pred[, arm] <- 1
        pred <- predict(mu_model, newdata = cbind(W_pred, X_test))$pred
        if (all(pred == 0)) {
          warning(paste("SuperLearner failed for mu at arm", arm, ". Using arm mean."))
          pred <- rep(mean(Y[in_train][W[in_train] == arm], na.rm = TRUE), n_test)
        }
        pred
      }),
      arm_levels
    )

    W_test_char <- as.character(W_test)
    mu_observed <- mapply(function(arm, j) mu_hat_per_arm[[arm]][j],
                          W_test_char, seq_along(W_test_char))

    pos <- setNames(
      lapply(all_pairs, function(pair) {
        arm_k <- pair[1]; arm_j <- pair[2]
        e_k   <- e_hat_test[, arm_k]
        e_j   <- e_hat_test[, arm_j]
        mu_k  <- mu_hat_per_arm[[arm_k]]
        mu_j  <- mu_hat_per_arm[[arm_j]]
        ind_k <- as.integer(W_test == arm_k)
        ind_j <- as.integer(W_test == arm_j)
        cate  <- mu_k - mu_j
        cate + (ind_k / e_k - ind_j / e_j) * (Y[in_test] - mu_observed)
      }),
      pair_names
    )

    e_per_arm <- setNames(lapply(arm_levels, function(arm) e_hat_test[, arm]), arm_levels)

    list(pos = pos, e_per_arm = e_per_arm, mu_hat_per_arm = mu_hat_per_arm, fold_pair = fold_pair)
  }, .options = furrr_options(seed = TRUE))

  po_matrices <- setNames(
    lapply(pair_names, function(pname) {
      collate_predictions(fold_list, fold_pairs, fold_indices,
                          lapply(cross_fits, function(cf) list(po = cf$pos[[pname]], fold_pair = cf$fold_pair)),
                          "po")
    }),
    pair_names
  )

  W.hat <- do.call(cbind, setNames(
    lapply(arm_levels, function(arm) {
      e_mat <- collate_predictions(fold_list, fold_pairs, fold_indices,
                                   lapply(cross_fits, function(cf) list(e_arm = cf$e_per_arm[[arm]], fold_pair = cf$fold_pair)),
                                   "e_arm")
      rowMeans(e_mat, na.rm = TRUE)
    }),
    arm_levels
  ))

  mu_per_arm <- setNames(
    lapply(arm_levels, function(arm) {
      mu_mat <- collate_predictions(fold_list, fold_pairs, fold_indices,
                                    lapply(cross_fits, function(cf) list(mu_arm = cf$mu_hat_per_arm[[arm]], fold_pair = cf$fold_pair)),
                                    "mu_arm")
      rowMeans(mu_mat, na.rm = TRUE)
    }),
    arm_levels
  )

  list(
    po_matrices = po_matrices,
    po          = setNames(lapply(pair_names, function(pn) rowMeans(po_matrices[[pn]], na.rm = TRUE)), pair_names),
    mu_per_arm  = mu_per_arm,
    W.hat       = W.hat
  )
}

# DR SuperLearner (multi-arm): stage_2_sl per pair
run_dr_superlearner_multiarm <- function(X, Y, W, nuisances, fold_indices, fold_list,
                                          arm_levels, all_pairs, sl_lib) {
  pair_names <- sapply(all_pairs, function(p) pair_name(p[1], p[2]))

  tau_list <- setNames(
    lapply(pair_names, function(pn) stage_2_sl(X, nuisances$po_matrices[[pn]], fold_indices, fold_list, sl_lib)),
    pair_names
  )

  list(
    tau         = tau_list,
    diagnostics = run_diagnostics_pairwise(X, Y, W, arm_levels, all_pairs, tau_list, nuisances)
  )
}

# Per-pair diagnostics: BLP (subset to pair, raw propensity) + independence test (all obs)
run_diagnostics_pairwise <- function(X, Y, W, arm_levels, all_pairs, tau_list, nuisances) {
  pair_names <- sapply(all_pairs, function(p) pair_name(p[1], p[2]))

  setNames(
    lapply(seq_along(all_pairs), function(p_idx) {
      pair  <- all_pairs[[p_idx]]
      arm_k <- pair[1]; arm_j <- pair[2]
      pname <- pair_names[p_idx]

      in_pair <- W %in% c(arm_k, arm_j)
      W_bin   <- as.integer(W[in_pair] == arm_k)
      e_k     <- nuisances$W.hat[in_pair, arm_k]
      Y0_j    <- nuisances$mu_per_arm[[arm_j]][in_pair]
      tau_sub <- tau_list[[pname]][in_pair]

      list(
        BLP_whole         = run_blp_whole(Y[in_pair], W_bin, e_k, Y0_j, tau_sub),
        independence_cate = run_independence_test_whole(X, tau_list[[pname]])
      )
    }),
    pair_names
  )
}
