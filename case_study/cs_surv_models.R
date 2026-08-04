##########
# Title: multi-arm CATE models for competing risks - factorial platform trial case study
# W must be a factor variable with K >= 2 levels; first level is the reference arm
# Returns results for ALL pairwise comparisons
##########
library(dplyr)
library(furrr)
library(grf)
library(pseudo)
library(SuperLearner)
library(ranger)
library(glmnet)
library(gam)

DEFAULT_SL_LIBRARY <- NULL

# ---- Utilities (self-contained) ----

collate_predictions <- function(
  fold_list,
  fold_pairs,
  fold_indices,
  reslist,
  target
) {
  lapply(fold_list, function(fold) {
    predictions <- rep(NA, length(fold_indices))
    for (j in seq_along(fold_pairs)) {
      if (fold %in% fold_pairs[[j]]) {
        predictions[fold_indices %in% fold_pairs[[j]]] <- reslist[[j]][[target]]
      }
    }
    predictions[fold_indices == fold] <- NA
    predictions
  }) %>%
    simplify2array()
}

stage_2_rf <- function(X, po, fold_indices, fold_list) {
  n_obs <- nrow(X)
  single <- is.vector(po)

  tau_results <- future_map(
    seq_along(fold_list),
    function(i) {
      fold <- fold_list[i]
      in_train <- fold_indices != fold
      in_fold <- !in_train
      forest <- if (single) {
        regression_forest(X[in_train, ], po[in_train])
      } else {
        regression_forest(X[in_train, ], po[in_train, fold])
      }
      list(
        fold = fold,
        predictions = predict(forest, newdata = X[in_fold, ])$predictions
      )
    },
    .options = furrr_options(seed = TRUE)
  )

  tau <- rep(NA_real_, n_obs)
  for (r in tau_results) {
    tau[fold_indices == r$fold] <- r$predictions
  }
  tau
}

pseudo_all <- function(Y, D, horizon) {
  D_int <- as.integer(D)
  Dc <- as.integer(D %in% c(1, 2))
  Y_sh <- ifelse(D == 2, horizon + 1, Y)
  D_sh <- as.integer(D == 1)

  ps_RMTL <- pseudoyl(Y, D_int, horizon)
  ps_RMSTc <- pseudomean(Y, Dc, horizon)
  ps_sh_RMST <- pseudomean(Y_sh, D_sh, horizon)

  list(
    ps_RMTL1 = ps_RMTL$pseudo$cause1,
    ps_RMTL2 = ps_RMTL$pseudo$cause2,
    ps_RMSTc = ps_RMSTc,
    ps_sh_RMST = ps_sh_RMST
  )
}

pseudo_crossfit <- function(
  Y,
  D,
  horizon,
  fold_indices,
  fold_list,
  pseudo_whole
) {
  n_obs <- length(Y)
  n_folds <- length(fold_list)
  D_int <- as.integer(D)
  Dc <- as.integer(D %in% c(1, 2))
  Y_sh <- ifelse(D == 2, horizon + 1, Y)
  D_sh <- as.integer(D == 1)

  pseudos <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold

    ps_RMTL <- pseudoyl(Y[in_train], D_int[in_train], horizon)
    ps_RMSTc <- pseudomean(Y[in_train], Dc[in_train], horizon)
    ps_sh_RMST <- pseudomean(Y_sh[in_train], D_sh[in_train], horizon)

    ps_RMTL$pseudo$cause1 <- ifelse(
      is.na(ps_RMTL$pseudo$cause1),
      pseudo_whole$ps_RMTL1[in_train],
      ps_RMTL$pseudo$cause1
    )
    ps_RMTL$pseudo$cause2 <- ifelse(
      is.na(ps_RMTL$pseudo$cause2),
      pseudo_whole$ps_RMTL2[in_train],
      ps_RMTL$pseudo$cause2
    )
    ps_RMSTc <- ifelse(
      is.na(ps_RMSTc),
      pseudo_whole$ps_RMSTc[in_train],
      ps_RMSTc
    )
    ps_sh_RMST <- ifelse(
      is.na(ps_sh_RMST),
      pseudo_whole$ps_sh_RMST[in_train],
      ps_sh_RMST
    )

    list(
      ps_RMTL1 = ps_RMTL$pseudo$cause1,
      ps_RMTL2 = ps_RMTL$pseudo$cause2,
      ps_RMSTc = ps_RMSTc,
      ps_sh_RMST = ps_sh_RMST,
      in_train = which(in_train),
      fold = i
    )
  })

  ps_RMTL1_mat <- matrix(NA_real_, nrow = n_obs, ncol = n_folds)
  ps_RMTL2_mat <- matrix(NA_real_, nrow = n_obs, ncol = n_folds)
  ps_RMSTc_mat <- matrix(NA_real_, nrow = n_obs, ncol = n_folds)
  ps_sh_RMST_mat <- matrix(NA_real_, nrow = n_obs, ncol = n_folds)

  for (result in pseudos) {
    i <- result$fold
    idx <- result$in_train
    ps_RMTL1_mat[idx, i] <- result$ps_RMTL1
    ps_RMTL2_mat[idx, i] <- result$ps_RMTL2
    ps_RMSTc_mat[idx, i] <- result$ps_RMSTc
    ps_sh_RMST_mat[idx, i] <- result$ps_sh_RMST
  }

  list(
    ps_RMTL1 = ps_RMTL1_mat,
    ps_RMTL2 = ps_RMTL2_mat,
    ps_RMSTc = ps_RMSTc_mat,
    ps_sh_RMST = ps_sh_RMST_mat
  )
}

pseudo_double_crossfit <- function(
  Y,
  D,
  horizon,
  fold_indices,
  fold_pairs,
  pseudo_whole
) {
  n_obs <- length(Y)
  D_int <- as.integer(D)
  Dc <- as.integer(D %in% c(1, 2))
  Y_sh <- ifelse(D == 2, horizon + 1, Y)
  D_sh <- as.integer(D == 1)

  pseudos <- future_map(seq_along(fold_pairs), function(i) {
    fold_pair <- fold_pairs[[i]]
    in_train <- !(fold_indices %in% fold_pair)

    ps_RMTL <- pseudoyl(Y[in_train], D_int[in_train], horizon)
    ps_RMSTc <- pseudomean(Y[in_train], Dc[in_train], horizon)
    ps_sh_RMST <- pseudomean(Y_sh[in_train], D_sh[in_train], horizon)

    ps_RMTL$pseudo$cause1 <- ifelse(
      is.na(ps_RMTL$pseudo$cause1),
      pseudo_whole$ps_RMTL1[in_train],
      ps_RMTL$pseudo$cause1
    )
    ps_RMTL$pseudo$cause2 <- ifelse(
      is.na(ps_RMTL$pseudo$cause2),
      pseudo_whole$ps_RMTL2[in_train],
      ps_RMTL$pseudo$cause2
    )
    ps_RMSTc <- ifelse(
      is.na(ps_RMSTc),
      pseudo_whole$ps_RMSTc[in_train],
      ps_RMSTc
    )
    ps_sh_RMST <- ifelse(
      is.na(ps_sh_RMST),
      pseudo_whole$ps_sh_RMST[in_train],
      ps_sh_RMST
    )

    list(
      ps_RMTL1 = ps_RMTL$pseudo$cause1,
      ps_RMTL2 = ps_RMTL$pseudo$cause2,
      ps_RMSTc = ps_RMSTc,
      ps_sh_RMST = ps_sh_RMST,
      fold_pair = fold_pair
    )
  })

  ps_RMTL1_mat <- matrix(NA_real_, nrow = n_obs, ncol = length(fold_pairs))
  ps_RMTL2_mat <- matrix(NA_real_, nrow = n_obs, ncol = length(fold_pairs))
  ps_RMSTc_mat <- matrix(NA_real_, nrow = n_obs, ncol = length(fold_pairs))
  ps_sh_RMST_mat <- matrix(NA_real_, nrow = n_obs, ncol = length(fold_pairs))

  for (i in seq_along(fold_pairs)) {
    in_train <- !(fold_indices %in% fold_pairs[[i]])
    ps_RMTL1_mat[in_train, i] <- pseudos[[i]]$ps_RMTL1
    ps_RMTL2_mat[in_train, i] <- pseudos[[i]]$ps_RMTL2
    ps_RMSTc_mat[in_train, i] <- pseudos[[i]]$ps_RMSTc
    ps_sh_RMST_mat[in_train, i] <- pseudos[[i]]$ps_sh_RMST
  }

  list(
    ps_RMTL1 = ps_RMTL1_mat,
    ps_RMTL2 = ps_RMTL2_mat,
    ps_RMSTc = ps_RMSTc_mat,
    ps_sh_RMST = ps_sh_RMST_mat
  )
}

# ---- Helpers ----

# Sign convention: pair_name(arm_k, arm_j) = tau(arm_k) - tau(arm_j)
# combn() pairs alphabetically, so "A1_vs_A2" = A1 - A2 (negate for active-vs-control)
pair_name <- function(arm_k, arm_j) paste0(arm_k, "_vs_", arm_j)

# Derive all-pairwise contrasts from multi_arm_causal_forest vs-reference predictions.
# macf_preds: n x (K-1) matrix (or plain vector when K=2 and drop=TRUE was used)
derive_all_pairwise <- function(macf_preds, arm_levels, all_pairs) {
  if (!is.matrix(macf_preds)) {
    macf_preds <- matrix(
      macf_preds,
      ncol = 1,
      dimnames = list(NULL, paste0(arm_levels[2], " - ", arm_levels[1]))
    )
  }

  ref_arm <- arm_levels[1]
  active_arms <- arm_levels[-1]

  vs_ref <- setNames(
    lapply(active_arms, function(arm) {
      macf_preds[, paste0(arm, " - ", ref_arm)]
    }),
    active_arms
  )

  setNames(
    lapply(all_pairs, function(pair) {
      arm_k <- pair[1]
      arm_j <- pair[2]
      if (arm_k == ref_arm) {
        -vs_ref[[arm_j]]
      } else if (arm_j == ref_arm) {
        vs_ref[[arm_k]]
      } else {
        vs_ref[[arm_k]] - vs_ref[[arm_j]]
      }
    }),
    sapply(all_pairs, function(p) pair_name(p[1], p[2]))
  )
}

# Adapted get_ipw: converts factor W to one-hot dummies before survival_forest
get_ipw_multiarm <- function(X, Y, D, W, horizon, censor) {
  D_event <- as.numeric(D == censor)

  W_dummies <- model.matrix(~ W - 1)
  sf_censor <- survival_forest(cbind(X, W_dummies), Y, D_event)

  censor_prob <- predict(
    sf_censor,
    failure.times = pmin(Y, horizon),
    prediction.times = "time"
  )$predictions

  observed <- (D != censor)
  epsilon <- 1e-3
  ipw <- 1 / pmax(censor_prob, epsilon)

  list(observed = observed, ipw = ipw)
}

# ---- Main orchestrator ----

all_cate_surv_models_multiarm <- function(
  data,
  Y,
  D,
  n_folds = 10,
  horizon = 30,
  sl_library = NULL
) {
  X <- as.matrix(data[, !(names(data) %in% c(Y, D, "W"))])
  Y <- data[[Y]]
  D <- data[[D]]
  W <- data$W

  if (!is.factor(W)) {
    stop("W must be a factor variable with named levels")
  }

  arm_levels <- levels(W)
  n_arms <- length(arm_levels)
  all_pairs <- combn(arm_levels, 2, simplify = FALSE)
  n_obs <- nrow(X)

  fold_indices <- sort(seq(n_obs) %% n_folds) + 1
  fold_list <- unique(fold_indices)
  fold_pairs <- utils::combn(fold_list, 2, simplify = FALSE)

  # Pre-compute propensity per fold (shared across all single-fold methods)
  # e_full[[i]]: n x K matrix; training rows use OOB preds, test-fold rows use out-of-fold preds
  message("Pre-computing propensity scores (probability_forest) per fold...")
  prop_per_fold <- future_map(
    seq_along(fold_list),
    function(i) {
      fold <- fold_list[i]
      in_train <- fold_indices != fold
      in_fold <- !in_train

      pf <- probability_forest(X[in_train, ], W[in_train])

      e_full <- matrix(NA_real_, nrow = n_obs, ncol = n_arms)
      colnames(e_full) <- arm_levels
      e_full[in_train, ] <- predict(pf)$predictions # OOB for train obs
      e_full[in_fold, ] <- predict(pf, newdata = X[in_fold, ])$predictions # out-of-fold for test obs

      e_full
    },
    .options = furrr_options(seed = TRUE)
  )

  results <- list()

  message("Multi-arm Causal Forest IPW approaches (RMST1, RMST2, RMSTc)...")
  results$ipw <- list()
  results$ipw$RMST1 <- cf_ipw_multiarm(
    X,
    Y,
    D,
    W,
    horizon,
    fold_indices,
    fold_list,
    prop_per_fold,
    arm_levels,
    all_pairs,
    event = 1
  )
  results$ipw$RMST2 <- cf_ipw_multiarm(
    X,
    Y,
    D,
    W,
    horizon,
    fold_indices,
    fold_list,
    prop_per_fold,
    arm_levels,
    all_pairs,
    event = 2
  )
  results$ipw$RMSTc <- cf_ipw_multiarm(
    X,
    Y,
    D,
    W,
    horizon,
    fold_indices,
    fold_list,
    prop_per_fold,
    arm_levels,
    all_pairs,
    event = "composite"
  )

  message(
    "Causal Survival Forest pairwise cs approaches (RMST1, RMST2, RMSTc)..."
  )
  results$csf_cs <- list()
  results$csf_cs$RMST1 <- csf_cs_multiarm(
    X,
    Y,
    D,
    W,
    horizon,
    fold_indices,
    fold_list,
    prop_per_fold,
    arm_levels,
    all_pairs,
    event = 1
  )
  results$csf_cs$RMST2 <- csf_cs_multiarm(
    X,
    Y,
    D,
    W,
    horizon,
    fold_indices,
    fold_list,
    prop_per_fold,
    arm_levels,
    all_pairs,
    event = 2
  )
  results$csf_cs$RMSTc <- csf_cs_multiarm(
    X,
    Y,
    D,
    W,
    horizon,
    fold_indices,
    fold_list,
    prop_per_fold,
    arm_levels,
    all_pairs,
    event = "composite"
  )

  message("Causal Survival Forest pairwise sh approaches (RMST1, RMST2)...")
  results$csf_sh <- list()
  results$csf_sh$RMST1 <- csf_sh_multiarm(
    X,
    Y,
    D,
    W,
    horizon,
    fold_indices,
    fold_list,
    prop_per_fold,
    arm_levels,
    all_pairs,
    event = 1
  )
  results$csf_sh$RMST2 <- csf_sh_multiarm(
    X,
    Y,
    D,
    W,
    horizon,
    fold_indices,
    fold_list,
    prop_per_fold,
    arm_levels,
    all_pairs,
    event = 2
  )

  message("Pseudo-value approaches...")
  pseudo_whole <- pseudo_all(Y, D, horizon)

  message("Pseudo-value multi_arm_causal_forest (RMTL1, RMTL2, RMSTc)...")
  pseudo_cv <- pseudo_crossfit(
    Y,
    D,
    horizon,
    fold_indices,
    fold_list,
    pseudo_whole
  )
  results$pseudo_cf <- list()
  results$pseudo_cf$RMTL1 <- pseudo_cf_multiarm(
    X,
    pseudo_cv$ps_RMTL1,
    W,
    fold_indices,
    fold_list,
    prop_per_fold,
    arm_levels,
    all_pairs
  )
  results$pseudo_cf$RMTL2 <- pseudo_cf_multiarm(
    X,
    pseudo_cv$ps_RMTL2,
    W,
    fold_indices,
    fold_list,
    prop_per_fold,
    arm_levels,
    all_pairs
  )
  results$pseudo_cf$RMSTc <- pseudo_cf_multiarm(
    X,
    pseudo_cv$ps_RMSTc,
    W,
    fold_indices,
    fold_list,
    prop_per_fold,
    arm_levels,
    all_pairs
  )

  message("Multi-arm pseudo DR learner (RMTL1, RMTL2, RMSTc)...")
  pseudo_dr_cv <- pseudo_double_crossfit(
    Y,
    D,
    horizon,
    fold_indices,
    fold_pairs,
    pseudo_whole
  )
  results$pseudo_dr <- list()
  nuisance_RMTL1 <- nuisance_pseudo_rf_multiarm(
    X,
    pseudo_dr_cv$ps_RMTL1,
    W,
    pseudo_whole$ps_RMTL1,
    fold_indices,
    fold_pairs,
    arm_levels,
    all_pairs
  )
  nuisance_RMTL2 <- nuisance_pseudo_rf_multiarm(
    X,
    pseudo_dr_cv$ps_RMTL2,
    W,
    pseudo_whole$ps_RMTL2,
    fold_indices,
    fold_pairs,
    arm_levels,
    all_pairs
  )
  nuisance_RMSTc <- nuisance_pseudo_rf_multiarm(
    X,
    pseudo_dr_cv$ps_RMSTc,
    W,
    pseudo_whole$ps_RMSTc,
    fold_indices,
    fold_pairs,
    arm_levels,
    all_pairs
  )

  results$pseudo_dr$RMTL1 <- pseudo_dr_rf_multiarm(
    X,
    nuisance_RMTL1,
    fold_indices,
    fold_list,
    all_pairs
  )
  results$pseudo_dr$RMTL2 <- pseudo_dr_rf_multiarm(
    X,
    nuisance_RMTL2,
    fold_indices,
    fold_list,
    all_pairs
  )
  results$pseudo_dr$RMSTc <- pseudo_dr_rf_multiarm(
    X,
    nuisance_RMSTc,
    fold_indices,
    fold_list,
    all_pairs
  )

  if (!is.null(sl_library)) {
    message(
      "SuperLearner T-learner multi-arm (standard pseudo-obs) (RMTL1, RMTL2, RMSTc)..."
    )
    results$sl_t_std <- list()
    results$sl_t_std$RMTL1 <- pseudo_sl_t_standard_multiarm(
      X,
      pseudo_cv$ps_RMTL1,
      W,
      fold_indices,
      fold_list,
      arm_levels,
      all_pairs,
      sl_library
    )
    results$sl_t_std$RMTL2 <- pseudo_sl_t_standard_multiarm(
      X,
      pseudo_cv$ps_RMTL2,
      W,
      fold_indices,
      fold_list,
      arm_levels,
      all_pairs,
      sl_library
    )
    results$sl_t_std$RMSTc <- pseudo_sl_t_standard_multiarm(
      X,
      pseudo_cv$ps_RMSTc,
      W,
      fold_indices,
      fold_list,
      arm_levels,
      all_pairs,
      sl_library
    )

    message(
      "SuperLearner DR-learner multi-arm (standard pseudo-obs) (RMTL1, RMTL2, RMSTc)..."
    )
    results$sl_dr_std <- list()
    sl_nuisance_RMTL1 <- nuisance_pseudo_sl_multiarm(
      X,
      pseudo_dr_cv$ps_RMTL1,
      W,
      pseudo_whole$ps_RMTL1,
      fold_indices,
      fold_pairs,
      arm_levels,
      all_pairs,
      sl_library
    )
    sl_nuisance_RMTL2 <- nuisance_pseudo_sl_multiarm(
      X,
      pseudo_dr_cv$ps_RMTL2,
      W,
      pseudo_whole$ps_RMTL2,
      fold_indices,
      fold_pairs,
      arm_levels,
      all_pairs,
      sl_library
    )
    sl_nuisance_RMSTc <- nuisance_pseudo_sl_multiarm(
      X,
      pseudo_dr_cv$ps_RMSTc,
      W,
      pseudo_whole$ps_RMSTc,
      fold_indices,
      fold_pairs,
      arm_levels,
      all_pairs,
      sl_library
    )
    results$sl_dr_std$RMTL1 <- pseudo_dr_sl_multiarm(
      X,
      sl_nuisance_RMTL1,
      fold_indices,
      fold_list,
      all_pairs,
      sl_library
    )
    results$sl_dr_std$RMTL2 <- pseudo_dr_sl_multiarm(
      X,
      sl_nuisance_RMTL2,
      fold_indices,
      fold_list,
      all_pairs,
      sl_library
    )
    results$sl_dr_std$RMSTc <- pseudo_dr_sl_multiarm(
      X,
      sl_nuisance_RMSTc,
      fold_indices,
      fold_list,
      all_pairs,
      sl_library
    )
  }

  results$fold_indices <- fold_indices
  results$arm_levels <- arm_levels
  results$all_pairs <- all_pairs

  return(results)
}

# ---- Method functions ----

# IPW with multi_arm_causal_forest
cf_ipw_multiarm <- function(
  X,
  Y,
  D,
  W,
  horizon,
  fold_indices,
  fold_list,
  prop_per_fold,
  arm_levels,
  all_pairs,
  event = 1
) {
  n_obs <- nrow(X)
  pair_names <- sapply(all_pairs, function(p) pair_name(p[1], p[2]))

  fold_results <- future_map(
    seq_along(fold_list),
    function(i) {
      fold <- fold_list[i]
      in_train <- fold_indices != fold
      in_fold <- !in_train

      weights_0 <- get_ipw_multiarm(
        X[in_train, ],
        Y[in_train],
        D[in_train],
        W[in_train],
        horizon,
        0
      )

      if (event == "composite") {
        total_observed <- weights_0$observed
        sample_weights <- weights_0$ipw[total_observed]
      } else {
        all_events <- unique(D[D != 0])
        competing <- setdiff(all_events, event)
        weights_ce <- get_ipw_multiarm(
          X[in_train, ],
          Y[in_train],
          D[in_train],
          W[in_train],
          horizon,
          competing
        )
        total_observed <- weights_0$observed & weights_ce$observed
        sample_weights <- (weights_0$ipw * weights_ce$ipw)[total_observed]
      }

      train_idx <- which(in_train)
      include <- train_idx[total_observed]
      W.hat_include <- prop_per_fold[[i]][include, , drop = FALSE]

      forest <- multi_arm_causal_forest(
        X[include, ],
        pmin(Y[include], horizon),
        W[include],
        W.hat = W.hat_include,
        sample.weights = sample_weights
      )

      pred <- predict(forest, newdata = X[in_fold, ])$predictions # n_fold x (K-1)
      list(fold = fold, pred = pred)
    },
    .options = furrr_options(seed = TRUE)
  )

  # Reconstruct n-vectors per pair
  result <- setNames(
    lapply(pair_names, function(pn) rep(NA_real_, n_obs)),
    pair_names
  )
  for (r in fold_results) {
    in_fold <- fold_indices == r$fold
    pairwise <- derive_all_pairwise(r$pred, arm_levels, all_pairs)
    for (pn in pair_names) {
      result[[pn]][in_fold] <- pairwise[[pn]]
    }
  }

  result
}

# Pairwise binary causal_survival_forest — cause-specific (CE treated as censoring)
csf_cs_multiarm <- function(
  X,
  Y,
  D,
  W,
  horizon,
  fold_indices,
  fold_list,
  prop_per_fold,
  arm_levels,
  all_pairs,
  event = 1
) {
  n_obs <- nrow(X)

  D_event <- if (event == "composite") {
    as.numeric(D %in% c(1, 2))
  } else {
    as.numeric(D == event)
  }

  pair_names <- sapply(all_pairs, function(p) pair_name(p[1], p[2]))

  setNames(
    lapply(seq_along(all_pairs), function(p_idx) {
      pair <- all_pairs[[p_idx]]
      arm_k <- pair[1]
      arm_j <- pair[2]

      fold_results <- future_map(
        seq_along(fold_list),
        function(i) {
          fold <- fold_list[i]
          in_train <- fold_indices != fold & W %in% c(arm_k, arm_j)
          in_fold <- fold_indices == fold

          W_bin <- as.integer(W[in_train] == arm_k)
          W.hat_kj <- prop_per_fold[[i]][in_train, arm_k]

          forest <- causal_survival_forest(
            X[in_train, ],
            Y[in_train],
            W_bin,
            D_event[in_train],
            W.hat = W.hat_kj,
            target = "RMST",
            horizon = horizon
          )

          pred <- predict(forest, newdata = X[in_fold, ])
          list(fold = fold, tau = pred$predictions)
        },
        .options = furrr_options(seed = TRUE)
      )

      tau <- rep(NA_real_, n_obs)
      for (r in fold_results) {
        tau[fold_indices == r$fold] <- r$tau
      }
      tau
    }),
    pair_names
  )
}

# Pairwise binary causal_survival_forest — shift approach (CE kept in risk set)
csf_sh_multiarm <- function(
  X,
  Y,
  D,
  W,
  horizon,
  fold_indices,
  fold_list,
  prop_per_fold,
  arm_levels,
  all_pairs,
  event = 1
) {
  n_obs <- nrow(X)

  D_sh <- as.numeric(D == event)
  Y_sh <- ifelse(!(D %in% c(event, 0)), horizon + 1, Y)

  pair_names <- sapply(all_pairs, function(p) pair_name(p[1], p[2]))

  setNames(
    lapply(seq_along(all_pairs), function(p_idx) {
      pair <- all_pairs[[p_idx]]
      arm_k <- pair[1]
      arm_j <- pair[2]

      fold_results <- future_map(
        seq_along(fold_list),
        function(i) {
          fold <- fold_list[i]
          in_train <- fold_indices != fold & W %in% c(arm_k, arm_j)

          include <- which(in_train)
          sample_weights <- NULL

          cens <- any(D[include] == 0 & Y[include] < horizon)
          if (cens) {
            weights_0 <- get_ipw_multiarm(
              X[include, ],
              Y[include],
              D[include],
              W[include],
              horizon,
              0
            )
            total_observed <- weights_0$observed
            sample_weights <- weights_0$ipw[total_observed]
            include <- include[total_observed]
          }

          W_bin <- as.integer(W[include] == arm_k)
          W.hat_kj <- prop_per_fold[[i]][include, arm_k]

          forest <- causal_survival_forest(
            X[include, ],
            Y_sh[include],
            W_bin,
            D_sh[include],
            W.hat = W.hat_kj,
            target = "RMST",
            horizon = horizon,
            sample.weights = sample_weights
          )

          pred <- predict(forest, newdata = X[fold_indices == fold, ])
          list(fold = fold, tau = pred$predictions)
        },
        .options = furrr_options(seed = TRUE)
      )

      tau <- rep(NA_real_, n_obs)
      for (r in fold_results) {
        tau[fold_indices == r$fold] <- r$tau
      }
      tau
    }),
    pair_names
  )
}

# Pseudo-value causal forest using multi_arm_causal_forest
# pseudo_mat: n x n_folds matrix (from pseudo_crossfit)
pseudo_cf_multiarm <- function(
  X,
  pseudo_mat,
  W,
  fold_indices,
  fold_list,
  prop_per_fold,
  arm_levels,
  all_pairs
) {
  n_obs <- nrow(X)
  pair_names <- sapply(all_pairs, function(p) pair_name(p[1], p[2]))

  fold_results <- future_map(
    seq_along(fold_list),
    function(i) {
      fold <- fold_list[i]
      in_train <- fold_indices != fold
      in_fold <- !in_train

      W.hat_train <- prop_per_fold[[i]][in_train, , drop = FALSE]

      forest <- multi_arm_causal_forest(
        X[in_train, ],
        pseudo_mat[in_train, fold],
        W[in_train],
        W.hat = W.hat_train
      )

      pred <- predict(forest, newdata = X[in_fold, ])$predictions # n_fold x (K-1)
      pairwise <- derive_all_pairwise(pred, arm_levels, all_pairs)
      list(fold = fold, pairwise = pairwise)
    },
    .options = furrr_options(seed = TRUE)
  )

  result <- setNames(
    lapply(pair_names, function(pn) rep(NA_real_, n_obs)),
    pair_names
  )
  for (r in fold_results) {
    in_fold <- fold_indices == r$fold
    for (pn in pair_names) {
      result[[pn]][in_fold] <- r$pairwise[[pn]]
    }
  }
  result
}

# Multi-arm DR learner nuisance estimation
# pseudo_mat: n x n_fold_pairs matrix (from pseudo_double_crossfit)
# Returns named list of po_matrices (n x n_folds each), one per pair
nuisance_pseudo_rf_multiarm <- function(
  X,
  pseudo_mat,
  W,
  pseudo_whole,
  fold_indices,
  fold_pairs,
  arm_levels,
  all_pairs
) {
  n_obs <- nrow(X)
  n_arms <- length(arm_levels)
  pair_names <- sapply(all_pairs, function(p) pair_name(p[1], p[2]))
  fold_list <- sort(unique(fold_indices))

  # Dummy matrix with arm-level column names (consistent for train and prediction)
  W_dummies <- model.matrix(~ W - 1)
  colnames(W_dummies) <- arm_levels # rename from "WA1" etc. to "A1" etc.

  cross_fits <- future_map(
    seq_along(fold_pairs),
    function(i) {
      fold_pair <- fold_pairs[[i]]
      in_train <- !(fold_indices %in% fold_pair)
      in_test <- !in_train
      n_test <- sum(in_test)

      X_train <- X[in_train, ]
      X_test <- X[in_test, ]
      W_train <- W[in_train]
      W_test <- W[in_test]
      W_dum_train <- W_dummies[in_train, , drop = FALSE]

      # Propensity: fit within the double-crossfit fold (correct for double crossfitting)
      pf <- probability_forest(X_train, W_train)
      e_hat_test <- predict(pf, newdata = X_test)$predictions # n_test x K
      # Ensure column names match arm_levels
      if (
        is.null(colnames(e_hat_test)) ||
          !all(colnames(e_hat_test) == arm_levels)
      ) {
        colnames(e_hat_test) <- arm_levels
      }

      # Outcome regression: pseudo ~ arm_dummies + X
      mu_model <- regression_forest(
        cbind(W_dum_train, X_train),
        pseudo_mat[in_train, i]
      )

      # Predict mu_k at each arm for test obs
      mu_hat_per_arm <- setNames(
        lapply(arm_levels, function(arm) {
          W_pred <- matrix(0, nrow = n_test, ncol = n_arms)
          colnames(W_pred) <- arm_levels
          W_pred[, arm] <- 1
          predict(mu_model, newdata = cbind(W_pred, X_test))$predictions
        }),
        arm_levels
      )

      # mu under each obs's observed arm
      W_test_char <- as.character(W_test)
      mu_observed <- mapply(
        function(arm, j) mu_hat_per_arm[[arm]][j],
        W_test_char,
        seq_along(W_test_char)
      )

      pseudo_obs <- pseudo_whole[in_test]

      # DR pseudo-outcome per pairwise comparison
      pos <- setNames(
        lapply(all_pairs, function(pair) {
          arm_k <- pair[1]
          arm_j <- pair[2]

          e_k <- e_hat_test[, arm_k]
          e_j <- e_hat_test[, arm_j]
          mu_k <- mu_hat_per_arm[[arm_k]]
          mu_j <- mu_hat_per_arm[[arm_j]]
          ind_k <- as.integer(W_test == arm_k)
          ind_j <- as.integer(W_test == arm_j)

          cate <- mu_k - mu_j
          po <- cate + (ind_k / e_k - ind_j / e_j) * (pseudo_obs - mu_observed)
          po
        }),
        pair_names
      )

      list(pos = pos, fold_pair = fold_pair)
    },
    .options = furrr_options(seed = TRUE)
  )

  # Collate into n x n_folds matrices per pair
  setNames(
    lapply(pair_names, function(pname) {
      cross_fits_pair <- lapply(cross_fits, function(cf) {
        list(po = cf$pos[[pname]], fold_pair = cf$fold_pair)
      })
      collate_predictions(
        fold_list,
        fold_pairs,
        fold_indices,
        cross_fits_pair,
        "po"
      )
    }),
    pair_names
  )
}

# DR learner stage 2: regression_forest on pseudo-outcomes, per pair
pseudo_dr_rf_multiarm <- function(
  X,
  po_matrices,
  fold_indices,
  fold_list,
  all_pairs
) {
  pair_names <- sapply(all_pairs, function(p) pair_name(p[1], p[2]))
  setNames(
    lapply(pair_names, function(pname) {
      stage_2_rf(X, po_matrices[[pname]], fold_indices, fold_list)
    }),
    pair_names
  )
}

# ---- SuperLearner helpers ----

# Stage 2 SuperLearner for DR-learner (po is an n x n_folds matrix)
stage_2_sl <- function(X, po, fold_indices, fold_list, sl_library = NULL) {
  n_obs <- nrow(X)

  tau_results <- future_map(
    seq_along(fold_list),
    function(i) {
      fold <- fold_list[i]
      in_train <- fold_indices != fold
      in_fold <- !in_train

      sl <- SuperLearner(
        Y = po[in_train, fold],
        X = as.data.frame(X[in_train, , drop = FALSE]),
        SL.library = sl_library,
        cvControl = list(V = 5)
      )
      tau_pred <- as.numeric(
        predict(sl, newdata = as.data.frame(X[in_fold, , drop = FALSE]))$pred
      )
      list(fold = fold, predictions = tau_pred)
    },
    .options = furrr_options(seed = TRUE)
  )

  tau <- rep(NA_real_, n_obs)
  for (result in tau_results) {
    tau[fold_indices == result$fold] <- result$predictions
  }
  tau
}

# ---- SuperLearner T-learner: standard cross-fit pseudo-obs, multi-arm ----
# pseudo_cv: n x n_folds matrix (from pseudo_crossfit)
# Fits one SL per arm on each training fold; returns named list of n-vectors (one per pair)
pseudo_sl_t_standard_multiarm <- function(
  X,
  pseudo_cv,
  W,
  fold_indices,
  fold_list,
  arm_levels,
  all_pairs,
  sl_library = NULL
) {
  n_obs <- nrow(X)
  pair_names <- sapply(all_pairs, function(p) pair_name(p[1], p[2]))

  fold_results <- future_map(
    seq_along(fold_list),
    function(i) {
      fold <- fold_list[i]
      in_train <- fold_indices != fold
      in_fold <- !in_train

      X_train <- as.data.frame(X[in_train, , drop = FALSE])
      X_test <- as.data.frame(X[in_fold, , drop = FALSE])
      W_train <- W[in_train]
      pseudo_train <- pseudo_cv[in_train, i]

      # One SL per arm, trained on that arm's observations only
      sl_per_arm <- setNames(
        lapply(arm_levels, function(arm) {
          in_arm <- W_train == arm
          SuperLearner(
            Y = pseudo_train[in_arm],
            X = X_train[in_arm, , drop = FALSE],
            SL.library = sl_library,
            cvControl = list(V = 5)
          )
        }),
        arm_levels
      )

      # mu_k(x) for all test observations under each arm
      mu_hat <- setNames(
        lapply(arm_levels, function(arm) {
          as.numeric(predict(sl_per_arm[[arm]], newdata = X_test)$pred)
        }),
        arm_levels
      )

      # Pairwise CATE = mu_k - mu_j
      pairwise <- setNames(
        lapply(all_pairs, function(pair) mu_hat[[pair[1]]] - mu_hat[[pair[2]]]),
        pair_names
      )

      list(fold = fold, pairwise = pairwise)
    },
    .options = furrr_options(seed = TRUE)
  )

  result <- setNames(
    lapply(pair_names, function(pn) rep(NA_real_, n_obs)),
    pair_names
  )
  for (r in fold_results) {
    in_fold <- fold_indices == r$fold
    for (pn in pair_names) {
      result[[pn]][in_fold] <- r$pairwise[[pn]]
    }
  }
  result
}

# ---- SuperLearner DR-learner nuisances: standard cross-fit pseudo-obs, multi-arm ----
# Mirrors nuisance_pseudo_rf_multiarm but replaces regression_forest with SuperLearner.
# Propensity uses probability_forest (multinomial; SL does not natively support K>2 classes).
# Returns a named list of po_matrices (n x n_fold_pairs each), one per pairwise contrast.
nuisance_pseudo_sl_multiarm <- function(
  X,
  pseudo_mat,
  W,
  pseudo_whole,
  fold_indices,
  fold_pairs,
  arm_levels,
  all_pairs,
  sl_library = NULL
) {
  n_obs <- nrow(X)
  n_arms <- length(arm_levels)
  pair_names <- sapply(all_pairs, function(p) pair_name(p[1], p[2]))
  fold_list <- sort(unique(fold_indices))

  W_dummies <- model.matrix(~ W - 1)
  colnames(W_dummies) <- arm_levels

  cross_fits <- future_map(
    seq_along(fold_pairs),
    function(i) {
      fold_pair <- fold_pairs[[i]]
      in_train <- !(fold_indices %in% fold_pair)
      in_test <- !in_train
      n_test <- sum(in_test)

      X_train <- as.data.frame(X[in_train, , drop = FALSE])
      X_test_df <- as.data.frame(X[in_test, , drop = FALSE])
      W_train <- W[in_train]
      W_test <- W[in_test]
      W_dum_train <- as.data.frame(W_dummies[in_train, , drop = FALSE])

      # Propensity: probability_forest (multinomial)
      pf <- probability_forest(X[in_train, ], W_train)
      e_hat_test <- predict(pf, newdata = X[in_test, ])$predictions
      if (
        is.null(colnames(e_hat_test)) ||
          !all(colnames(e_hat_test) == arm_levels)
      ) {
        colnames(e_hat_test) <- arm_levels
      }

      # Outcome regression with SL: pseudo ~ arm_dummies + X
      sl_mu <- SuperLearner(
        Y = pseudo_mat[in_train, i],
        X = cbind(W_dum_train, X_train),
        SL.library = sl_library,
        cvControl = list(V = 5)
      )

      # mu_k(x) for each arm at test observations
      mu_hat_per_arm <- setNames(
        lapply(arm_levels, function(arm) {
          W_pred <- as.data.frame(matrix(
            0,
            nrow = n_test,
            ncol = n_arms,
            dimnames = list(NULL, arm_levels)
          ))
          W_pred[, arm] <- 1
          as.numeric(predict(sl_mu, newdata = cbind(W_pred, X_test_df))$pred)
        }),
        arm_levels
      )

      W_test_char <- as.character(W_test)
      mu_observed <- mapply(
        function(arm, j) mu_hat_per_arm[[arm]][j],
        W_test_char,
        seq_along(W_test_char)
      )
      pseudo_obs <- pseudo_whole[in_test]

      # DR pseudo-outcome per pairwise comparison
      pos <- setNames(
        lapply(all_pairs, function(pair) {
          arm_k <- pair[1]
          arm_j <- pair[2]
          e_k <- e_hat_test[, arm_k]
          e_j <- e_hat_test[, arm_j]
          mu_k <- mu_hat_per_arm[[arm_k]]
          mu_j <- mu_hat_per_arm[[arm_j]]
          ind_k <- as.integer(W_test == arm_k)
          ind_j <- as.integer(W_test == arm_j)
          cate <- mu_k - mu_j
          po <- cate + (ind_k / e_k - ind_j / e_j) * (pseudo_obs - mu_observed)
          po
        }),
        pair_names
      )

      list(pos = pos, fold_pair = fold_pair)
    },
    .options = furrr_options(seed = TRUE)
  )

  setNames(
    lapply(pair_names, function(pname) {
      cross_fits_pair <- lapply(cross_fits, function(cf) {
        list(po = cf$pos[[pname]], fold_pair = cf$fold_pair)
      })
      collate_predictions(
        fold_list,
        fold_pairs,
        fold_indices,
        cross_fits_pair,
        "po"
      )
    }),
    pair_names
  )
}

# DR learner stage 2: SuperLearner on pseudo-outcomes, per pair
pseudo_dr_sl_multiarm <- function(
  X,
  po_matrices,
  fold_indices,
  fold_list,
  all_pairs,
  sl_library = NULL
) {
  pair_names <- sapply(all_pairs, function(p) pair_name(p[1], p[2]))
  setNames(
    lapply(pair_names, function(pname) {
      stage_2_sl(X, po_matrices[[pname]], fold_indices, fold_list, sl_library)
    }),
    pair_names
  )
}
