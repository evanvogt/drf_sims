##########
# Title: running CATE models - Competing risks
##########
library(dplyr)
library(furrr)
library(grf)
library(pseudo)
library(SuperLearner)
library(ranger)
library(glmnet)
library(gam)

DEFAULT_SL_LIBRARY <- c("SL.glm", "SL.glmnet", "SL.ranger", "SL.gam")

all_cate_surv_models <- function(
  data,
  n_folds = 10,
  horizon = 30,
  sl_library = DEFAULT_SL_LIBRARY
) {
  # data formatting
  X <- as.matrix(data[, !(names(data) %in% c("Y", "D", "W"))])
  Y <- data$Y
  D <- data$D
  W <- data$W
  n_obs <- nrow(X)

  # Fold creation. Only the "scf" arms and the split-pseudo T-learner still split;
  # the whole_oob arms take grf's own OOB predictions and ignore these.
  fold_indices <- sort(seq(n_obs) %% n_folds) + 1
  fold_list <- unique(fold_indices)

  # all results container
  results <- list()

  message("Causal Forest IPW approaches (RMST1, RMST2, RMSTc)...")
  results$ipw <- list()
  results$ipw$RMST1 <- cf_ipw(X, Y, D, W, horizon, event = 1)
  results$ipw$RMST2 <- cf_ipw(X, Y, D, W, horizon, event = 2)
  results$ipw$RMSTc <- cf_ipw(X, Y, D, W, horizon, event = "composite")

  message("Causal Survival Forest cs approaches (RMST1, RMST2, RMSTc)...")
  results$csf_cs <- list()
  results$csf_cs$RMST1 <- csf_cs(X, Y, D, W, horizon, event = 1)
  results$csf_cs$RMST2 <- csf_cs(X, Y, D, W, horizon, event = 2)
  results$csf_cs$RMSTc <- csf_cs(X, Y, D, W, horizon, event = "composite")

  message("Causal Survival Forest sh approaches (RMST1, RMST2)...")
  results$csf_sh <- list()
  results$csf_sh$RMST1 <- csf_sh(X, Y, D, W, horizon, event = 1)
  results$csf_sh$RMST2 <- csf_sh(X, Y, D, W, horizon, event = 2)

  # ---- pseudo-value arms ----------------------------------------------------
  # Two factors, crossed as far as they can be (see README):
  #   pseudo-values  "whole" (one AJ/KM fit on all n) vs "cvps" (leave-one-fold-out)
  #   fitting        "oob" (whole sample, grf-internal) vs "scf" (single crossfit)
  # cvps + oob is not a cell: the crossfit matrix is NA on each fold's own rows,
  # so it only exists inside a fold loop. whole_scf is the control that separates
  # the pseudo-value factor from the fitting factor.
  estimands <- c("RMTL1", "RMTL2", "RMSTc")
  by_estimand <- function(f) setNames(lapply(estimands, f), estimands)

  message("Pseudo-value approaches...")
  pseudo_whole <- pseudo_all(Y, D, horizon)
  pseudo_cv <- pseudo_crossfit(
    Y,
    D,
    horizon,
    fold_indices,
    fold_list,
    pseudo_whole
  )
  ps_whole <- function(e) pseudo_whole[[paste0("ps_", e)]]
  ps_cv <- function(e) pseudo_cv[[paste0("ps_", e)]]

  message("Pseudo-value Causal Forest (whole_oob, whole_scf, cvps_scf)...")
  results$pseudo_cf_whole_oob <- by_estimand(function(e) {
    pseudo_cf_whole_oob(X, ps_whole(e), W)
  })
  results$pseudo_cf_whole_scf <- by_estimand(function(e) {
    pseudo_cf_scf(X, ps_whole(e), W, fold_indices, fold_list)
  })
  results$pseudo_cf_cvps_scf <- by_estimand(function(e) {
    pseudo_cf_scf(X, ps_cv(e), W, fold_indices, fold_list)
  })

  message("Pseudo-value DR Learner nuisances (whole_oob, whole_scf, cvps_scf)...")
  nuis_dr <- list(
    whole_oob = by_estimand(function(e) {
      nuisance_pseudo_rf_oob(X, ps_whole(e), W)
    }),
    whole_scf = by_estimand(function(e) {
      nuisance_pseudo_rf_scf(
        X,
        ps_whole(e),
        W,
        ps_whole(e),
        fold_indices,
        fold_list
      )
    }),
    cvps_scf = by_estimand(function(e) {
      nuisance_pseudo_rf_scf(
        X,
        ps_cv(e),
        W,
        ps_whole(e),
        fold_indices,
        fold_list
      )
    })
  )

  message("DR second stage regression (RMTL1, RMTL2, RMSTc)...")
  results$pseudo_dr_whole_oob <- by_estimand(function(e) {
    stage2_whole_rf(X, nuis_dr$whole_oob[[e]]$po)$tau
  })
  results$pseudo_dr_whole_scf <- by_estimand(function(e) {
    stage_2_rf_scf(X, nuis_dr$whole_scf[[e]]$po, fold_indices, fold_list)
  })
  results$pseudo_dr_cvps_scf <- by_estimand(function(e) {
    stage_2_rf_scf(X, nuis_dr$cvps_scf[[e]]$po, fold_indices, fold_list)
  })

  message("SuperLearner T-learner (whole and cvps pseudo-obs)...")
  results$sl_t_whole <- by_estimand(function(e) {
    pseudo_sl_t_standard(X, ps_whole(e), W, fold_indices, fold_list, sl_library)
  })
  results$sl_t_cvps <- by_estimand(function(e) {
    pseudo_sl_t_standard(X, ps_cv(e), W, fold_indices, fold_list, sl_library)
  })

  message("SuperLearner T-learner (split pseudo-obs)...")
  results$sl_t_split <- pseudo_sl_t_split(
    X,
    Y,
    D,
    W,
    horizon,
    fold_indices,
    fold_list,
    sl_library
  )

  message("SuperLearner DR-learner (whole and cvps pseudo-obs)...")
  nuis_sl <- list(
    whole = by_estimand(function(e) {
      nuisance_pseudo_sl(
        X,
        ps_whole(e),
        W,
        ps_whole(e),
        fold_indices,
        fold_list,
        sl_library
      )
    }),
    cvps = by_estimand(function(e) {
      nuisance_pseudo_sl(
        X,
        ps_cv(e),
        W,
        ps_whole(e),
        fold_indices,
        fold_list,
        sl_library
      )
    })
  )
  results$sl_dr_whole <- by_estimand(function(e) {
    pseudo_dr_sl(X, nuis_sl$whole[[e]]$po, fold_indices, fold_list, sl_library)
  })
  results$sl_dr_cvps <- by_estimand(function(e) {
    pseudo_dr_sl(X, nuis_sl$cvps[[e]]$po, fold_indices, fold_list, sl_library)
  })

  # SuperLearner DR-learner (split pseudo-obs) - DISABLED for now, see
  # README.md "Known issues": compute_split_pseudoyl()/compute_split_pseudomean()
  # (below) hit an unpatched NaN in pseudo::pseudoyl()'s internal ci.omit()
  # whenever a validation Y exceeds the max Y in that fold's KM set, and
  # unlike pseudo_crossfit() there is no NA
  # fallback here, so the NA reaches stage_2_sl_vec()'s SuperLearner() call
  # and aborts the run. Diagnosed and reproduced in
  # surv_dr_split_na_diagnose.R. nuisance_pseudo_sl_split()/stage_2_sl_vec()
  # are left defined below for that diagnostic script and for the fix.
  #
  # message("SuperLearner DR-learner (split pseudo-obs)...")
  # sl_split_nuisances <- nuisance_pseudo_sl_split(
  #   X,
  #   Y,
  #   D,
  #   W,
  #   horizon,
  #   fold_indices,
  #   fold_list,
  #   sl_library
  # )
  # results$sl_dr_split <- list()
  # results$sl_dr_split$RMTL1 <- stage_2_sl_vec(
  #   X,
  #   sl_split_nuisances$po_RMTL1,
  #   fold_indices,
  #   fold_list,
  #   sl_library
  # )
  # results$sl_dr_split$RMTL2 <- stage_2_sl_vec(
  #   X,
  #   sl_split_nuisances$po_RMTL2,
  #   fold_indices,
  #   fold_list,
  #   sl_library
  # )
  # results$sl_dr_split$RMSTc <- stage_2_sl_vec(
  #   X,
  #   sl_split_nuisances$po_RMSTc,
  #   fold_indices,
  #   fold_list,
  #   sl_library
  # )

  # add extra bits to the results
  results$pseudos <- list()
  results$pseudos$cf_cv <- pseudo_cv
  results$pseudos$whole <- pseudo_whole

  # Nuisances are now length-n vectors rather than n x C(V,2) matrices, so there
  # is nothing left to aggregate with rowMeans - the same simplification
  # R/cate_models.R made when it dropped its aggregate_nuisances axis. Nested by
  # arm, then estimand.
  results$nuisances <- list(rf = nuis_dr, sl = nuis_sl)

  results$fold_indices <- fold_indices

  return(results)
}

# Helper Functions
# get IPW weights for specific event (censoring or the competing event)
get_ipw <- function(X, Y, D, W, horizon, censor) {
  # select censoring event
  D_event <- as.numeric(D == censor)

  # fit survival forest for time to censoring event
  sf_censor <- survival_forest(cbind(X, W), Y, D_event)

  censor_prob <- predict(
    sf_censor,
    failure.times = pmin(Y, horizon),
    prediction.times = "time"
  )$predictions

  # get weights for non censored (clipped)
  observed <- (D != censor)
  epsilon <- 1e-3
  ipw <- 1 / pmax(censor_prob, epsilon)

  list(observed = observed, ipw = ipw)
}
# CF using IPW weights for censoring and CEs - see GRF tutorial
#
# Whole-sample, grf-internal crossfitting ("cf_default"), matching
# R/cate_models.R::run_causal_forest. The forest is fit only on `include` -
# observations not censored and free of the competing event - so OOB predictions
# exist for those rows only; the excluded rows never entered the forest at all,
# so a plain newdata prediction for them is honest. Same pattern as
# crossfitting/cf_models.R::nuisance_oob_rf.
cf_ipw <- function(X, Y, D, W, horizon, event = 1) {
  n_obs <- nrow(X)

  # get ipw weights. get_ipw()'s predict() call passes no newdata, so the
  # censoring probabilities are already grf OOB predictions.
  weights_0 <- get_ipw(X, Y, D, W, horizon, 0) # all 1's when there is no censoring

  if (event == "composite") {
    total_observed <- weights_0$observed
    sample_weights <- weights_0$ipw[total_observed]
  } else {
    # Identify event of interest
    all_events <- unique(D[D != 0])
    competing <- setdiff(all_events, event)

    weights_competing <- get_ipw(X, Y, D, W, horizon, competing)

    total_observed <- weights_0$observed & weights_competing$observed
    sample_weights <- weights_0$ipw * weights_competing$ipw
    sample_weights <- sample_weights[total_observed]
  }

  include <- which(total_observed)

  # causal forest (not survival because we've sorted out weights - see grf tutorial)
  # Y truncated at horizon so the outcome is min(T1*, horizon), matching the RMST estimand
  forest <- causal_forest(
    X[include, ],
    pmin(Y[include], horizon),
    W[include],
    sample.weights = sample_weights
  )

  tau_RMST <- rep(NA_real_, n_obs)
  tau_RMST[include] <- predict(forest)$predictions # OOB for the rows it was fit on
  if (length(include) < n_obs) {
    tau_RMST[-include] <- predict(
      forest,
      newdata = X[-include, , drop = FALSE]
    )$predictions
  }

  return(tau_RMST)
}
# CSF - treating CEs as censoring events
# Whole sample; causal_survival_forest cross-fits its own nuisances internally
# and predict() with no newdata returns OOB predictions.
csf_cs <- function(X, Y, D, W, horizon, event = 1) {
  if (event == "composite") {
    # Modify event to include 1 and 2
    D_event <- as.numeric(D %in% c(1, 2))
  } else {
    # Modify event def to treat competing event as censoring
    D_event <- as.numeric(D == event)
  }

  forest <- causal_survival_forest(
    X,
    Y,
    W,
    D_event,
    target = "RMST",
    horizon = horizon
  )

  return(predict(forest)$predictions)
}
# CSF - keep CEs in the risk set
# Whole sample, as csf_cs. When censoring is present the forest is fit on the
# uncensored subset only, so the excluded rows get a newdata prediction - see
# cf_ipw above for why that is honest.
csf_sh <- function(X, Y, D, W, horizon, event = 1) {
  n_obs <- nrow(X)

  # move competing events after horizon (keep them in the risk set)
  D_sh <- as.numeric(D == event)
  Y_sh <- ifelse(!(D %in% c(event, 0)), horizon + 1, Y)

  include <- seq_len(n_obs)
  sample_weights <- NULL
  # if there is censoring for event horizon, account for this
  cens <- any(D == 0 & Y < horizon)
  if (cens) {
    weights_0 <- get_ipw(X, Y, D, W, horizon, 0)
    include <- which(weights_0$observed)
    sample_weights <- weights_0$ipw[weights_0$observed]
  }

  forest <- causal_survival_forest(
    X[include, ],
    Y_sh[include],
    W[include],
    D_sh[include],
    target = "RMST",
    horizon = horizon,
    sample.weights = sample_weights
  )

  tau_RMST <- rep(NA_real_, n_obs)
  tau_RMST[include] <- predict(forest)$predictions
  if (length(include) < n_obs) {
    tau_RMST[-include] <- predict(
      forest,
      newdata = X[-include, , drop = FALSE]
    )$predictions
  }

  return(tau_RMST)
}
# pseudovalue aproaches
#' Jackknife pseudo-values from one Aalen-Johansen / Kaplan-Meier fit on all n rows
#'
#' NOTE: a fourth estimand, `ps_sh_RMST` (subdistribution / Fine-Gray RMST, on
#' `Y_sh = ifelse(D == 2, horizon + 1, Y)` and `D_sh = as.integer(D == 1)`), used
#' to be computed here and in every crossfit variant below. It was the pseudo-value
#' counterpart to the `csf_sh` framework, but no pseudo-value arm on the
#' subdistribution scale was ever built, so across the whole of this file's history
#' it was computed and stored and consumed by no estimator. Removed rather than
#' carried; reinstate it here if a subdistribution pseudo-value arm is added.
pseudo_all <- function(Y, D, horizon) {
  # reformatting for different estimands
  D_int <- as.integer(D)
  Dc <- as.integer(D %in% c(1, 2))

  ps_RMTL <- pseudoyl(Y, D_int, horizon)
  ps_RMSTc <- pseudomean(Y, Dc, horizon)

  list(
    ps_RMTL1 = ps_RMTL$pseudo$cause1,
    ps_RMTL2 = ps_RMTL$pseudo$cause2,
    ps_RMSTc = ps_RMSTc
  )
}
#' Leave-one-fold-out pseudo-values ("cvps")
#'
#' Column k holds pseudo-values from an AJ/KM fit that never saw fold k, and is
#' NA on fold k's own rows - a jackknife pseudo-value for observation j requires
#' j to be in the sample being decomposed, so there is no honest out-of-fold value
#' for the held-out rows. That is why the "cvps" arms only exist inside a fold
#' loop, and why the DR arms keep whole-sample pseudo-values in their correction
#' term (see the README).
#'
#' `n_na_fallback` records how often pseudoyl()/pseudomean() returned NA and the
#' whole-sample value was substituted. That substitution leaks the held-out fold
#' into the "cvps" arms, so the count qualifies the whole-vs-crossfit comparison
#' rather than being a diagnostic to ignore. Per README "Known issues" it fires on
#' the max-time observation of each fold, so expect it to be O(V).
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

  # reformatting for different estimands
  D_int <- as.integer(D)
  Dc <- as.integer(D %in% c(1, 2))

  pseudos <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold

    ps_RMTL <- pseudoyl(Y[in_train], D_int[in_train], horizon)
    ps_RMSTc <- pseudomean(Y[in_train], Dc[in_train], horizon)

    ps_RMTL1 <- ps_RMTL$pseudo$cause1
    ps_RMTL2 <- ps_RMTL$pseudo$cause2

    # count before substituting, so the leakage it introduces stays measurable
    n_na <- c(
      RMTL1 = sum(is.na(ps_RMTL1)),
      RMTL2 = sum(is.na(ps_RMTL2)),
      RMSTc = sum(is.na(ps_RMSTc))
    )

    ps_RMTL1 <- ifelse(
      is.na(ps_RMTL1),
      pseudo_whole$ps_RMTL1[in_train],
      ps_RMTL1
    )
    ps_RMTL2 <- ifelse(
      is.na(ps_RMTL2),
      pseudo_whole$ps_RMTL2[in_train],
      ps_RMTL2
    )
    ps_RMSTc <- ifelse(
      is.na(ps_RMSTc),
      pseudo_whole$ps_RMSTc[in_train],
      ps_RMSTc
    )

    list(
      ps_RMTL1 = ps_RMTL1,
      ps_RMTL2 = ps_RMTL2,
      ps_RMSTc = ps_RMSTc,
      n_na = n_na,
      in_train = which(in_train),
      fold = i
    )
  })

  # Empty matrices for pseudos
  ps_RMTL1_mat <- matrix(NA_real_, nrow = n_obs, ncol = n_folds)
  ps_RMTL2_mat <- matrix(NA_real_, nrow = n_obs, ncol = n_folds)
  ps_RMSTc_mat <- matrix(NA_real_, nrow = n_obs, ncol = n_folds)

  # Fill matrices
  for (result in pseudos) {
    i <- result$fold
    idx <- result$in_train
    ps_RMTL1_mat[idx, i] <- result$ps_RMTL1
    ps_RMTL2_mat[idx, i] <- result$ps_RMTL2
    ps_RMSTc_mat[idx, i] <- result$ps_RMSTc
  }
  list(
    ps_RMTL1 = ps_RMTL1_mat,
    ps_RMTL2 = ps_RMTL2_mat,
    ps_RMSTc = ps_RMSTc_mat,
    n_na_fallback = colSums(do.call(rbind, lapply(pseudos, `[[`, "n_na")))
  )
}
# causal forest using pseudo values - whole sample, grf's own internal
# crossfitting ("cf_default"), matching R/cate_models.R::run_causal_forest.
# `pseudo` is the whole-sample pseudo-value vector.
pseudo_cf_whole_oob <- function(X, pseudo, W) {
  forest <- causal_forest(X, pseudo, W)
  predict(forest)$predictions
}

# causal forest using pseudo values - single leave-one-fold-out ("scf").
#
# `pseudo` is either the whole-sample vector (the whole_scf control arm) or the
# n x V crossfit matrix from pseudo_crossfit (the cvps_scf arm). Branching on
# is.matrix() here is what lets one function serve both arms, so that the two
# differ in the pseudo-values alone and in nothing else.
pseudo_cf_scf <- function(X, pseudo, W, fold_indices, fold_list) {
  n_obs <- nrow(X)
  cvps <- is.matrix(pseudo)

  tau_result <- future_map(
    seq_along(fold_list),
    function(i) {
      fold <- fold_list[i]
      in_train <- fold_indices != fold
      in_fold <- !in_train

      pseudo_train <- if (cvps) pseudo[in_train, fold] else pseudo[in_train]

      forest <- causal_forest(
        X[in_train, ],
        pseudo_train,
        W[in_train]
      )

      pred <- predict(forest, newdata = X[in_fold, ])

      list(fold = fold, tau = pred$predictions)
    },
    .options = furrr_options(seed = TRUE)
  )

  tau <- rep(NA_real_, n_obs)
  for (result in tau_result) {
    in_fold <- fold_indices == result$fold
    tau[in_fold] <- result$tau
  }
  return(tau)
}
# DR learner nuisance functions
#
# Whole-sample OOB, S-learner ("oob_oob_s") - the pseudo-value counterpart of
# R/cate_models.R::nuisance_rf, and the production-parity arm. No sample
# splitting: one S-learner forest on cbind(W, X) supplies both counterfactuals
# via oob_predict_counterfactual, and two more whole-sample forests supply W.hat
# and pseudo.hat.cf, all taken out-of-bag.
#
# `pseudo` is the whole-sample pseudo-value vector, and it is also the outcome in
# the DR correction term - there is no separate `pseudo_whole` argument here
# because with no split the two coincide.
nuisance_pseudo_rf_oob <- function(X, pseudo, W) {
  forest <- regression_forest(cbind(W = W, X), pseudo)
  pseudo0.hat <- oob_predict_counterfactual(forest, cbind(W = 0, X))
  pseudo1.hat <- oob_predict_counterfactual(forest, cbind(W = 1, X))

  W.hat <- trim_ps(predict(regression_forest(X, W))$predictions)
  pseudo.hat.cf <- predict(regression_forest(X, pseudo))$predictions

  pseudo.hat <- W * pseudo1.hat + (1 - W) * pseudo0.hat
  po <- dr_pseudo(pseudo, W, pseudo1.hat, pseudo0.hat, W.hat)

  list(
    po = po,
    pseudo.hat = pseudo.hat,
    pseudo0.hat = pseudo0.hat,
    pseudo.hat.cf = pseudo.hat.cf,
    W.hat = W.hat
  )
}

# DR learner nuisances - single leave-one-fold-out ("scf").
#
# `pseudo` is either the whole-sample vector (whole_scf) or the n x V crossfit
# matrix (cvps_scf); the is.matrix() branch is the only difference between those
# two arms.
#
# `pseudo_whole` stays the outcome in the correction term for BOTH arms. This is
# deliberate, not an oversight: pseudo_crossfit has no value for the held-out
# rows (a jackknife pseudo-value for observation j needs j in the decomposed
# sample), so the factor these arms vary is the pseudo-values used to TRAIN the
# nuisance regressions, not the ones entering the DR correction. See the README.
nuisance_pseudo_rf_scf <- function(
  X,
  pseudo,
  W,
  pseudo_whole,
  fold_indices,
  fold_list
) {
  cvps <- is.matrix(pseudo)

  cross_fits <- future_map(
    seq_along(fold_list),
    function(i) {
      fold <- fold_list[i]
      in_train <- fold_indices != fold
      in_test <- !in_train

      pseudo_train <- if (cvps) pseudo[in_train, fold] else pseudo[in_train]

      ps.hat.model <- regression_forest(
        cbind(W[in_train], X[in_train, ]),
        pseudo_train
      )
      ps.hat.cf.model <- regression_forest(X[in_train, ], pseudo_train)
      W.hat.model <- regression_forest(X[in_train, ], W[in_train])

      X_test <- X[in_test, ]

      pseudo0.hat <- predict(
        ps.hat.model,
        newdata = cbind(W = 0, X_test)
      )$predictions
      pseudo1.hat <- predict(
        ps.hat.model,
        newdata = cbind(W = 1, X_test)
      )$predictions
      pseudo.hat.cf <- predict(ps.hat.cf.model, newdata = X_test)$predictions
      W.hat <- trim_ps(predict(W.hat.model, newdata = X_test)$predictions)

      # DR learner pseudo outcome (not the same as the pseudo values we already have)
      W_test <- W[in_test]
      pseudo.hat <- W_test * pseudo1.hat + (1 - W_test) * pseudo0.hat
      po <- dr_pseudo(
        pseudo_whole[in_test],
        W_test,
        pseudo1.hat,
        pseudo0.hat,
        W.hat
      )

      list(
        fold = fold,
        po = po,
        pseudo.hat = pseudo.hat,
        pseudo0.hat = pseudo0.hat,
        pseudo.hat.cf = pseudo.hat.cf,
        W.hat = W.hat
      )
    },
    .options = furrr_options(seed = TRUE)
  )

  scatter_folds(
    cross_fits,
    fold_indices,
    c("po", "pseudo.hat", "pseudo0.hat", "pseudo.hat.cf", "W.hat")
  )
}

# Leave-one-fold-out second stage on a po VECTOR. R/cate_models.R no longer has a
# fold-wise RF stage 2 to borrow (it moved to whole-sample OOB, stage2_whole_rf),
# so the scf arms keep this local one. The whole_oob arm uses stage2_whole_rf.
stage_2_rf_scf <- function(X, po, fold_indices, fold_list) {
  n_obs <- nrow(X)

  tau_results <- future_map(
    seq_along(fold_list),
    function(i) {
      fold <- fold_list[i]
      in_train <- fold_indices != fold
      in_fold <- !in_train

      forest <- regression_forest(X[in_train, ], po[in_train])

      tau_pred <- predict(forest, newdata = X[in_fold, ])$predictions

      list(fold = fold, predictions = tau_pred)
    },
    .options = furrr_options(seed = TRUE)
  )

  # Reconstruct tau vector
  tau <- rep(NA_real_, n_obs)
  for (result in tau_results) {
    tau[fold_indices == result$fold] <- result$predictions
  }
  return(tau)
}

# Split pseudo-observation helpers (Cwiling et al. 2025, Eq. 3)
# Key identity: split pseudo-obs for obs i added to KM set of size n_km equals
# the last element of pseudomean(c(Y_km, Y[i]), c(D_km, D[i]), horizon),
# because removing obs i from the combined dataset recovers D_{n1}.
compute_split_pseudomean <- function(Y_km, D_km, Y_val, D_val, horizon) {
  n_km <- length(Y_km)
  vapply(
    seq_along(Y_val),
    function(j) {
      pseudomean(c(Y_km, Y_val[j]), c(D_km, D_val[j]), horizon)[n_km + 1]
    },
    numeric(1)
  )
}

compute_split_pseudoyl <- function(Y_km, D_km, Y_val, D_val, horizon) {
  n_km <- length(Y_km)
  cause1 <- vapply(
    seq_along(Y_val),
    function(j) {
      pseudoyl(
        c(Y_km, Y_val[j]),
        as.integer(c(D_km, D_val[j])),
        horizon
      )$pseudo$cause1[n_km + 1]
    },
    numeric(1)
  )
  cause2 <- vapply(
    seq_along(Y_val),
    function(j) {
      pseudoyl(
        c(Y_km, Y_val[j]),
        as.integer(c(D_km, D_val[j])),
        horizon
      )$pseudo$cause2[n_km + 1]
    },
    numeric(1)
  )
  list(cause1 = cause1, cause2 = cause2)
}

# SuperLearner T-learner, single leave-one-fold-out ("scf_scf").
#
# `pseudo` is either the whole-sample vector (sl_t_whole) or the n x V crossfit
# matrix from pseudo_crossfit (sl_t_cvps). SuperLearner has no OOB analogue - see
# crossfitting/README.md, which drops the OOB arms for the SL comparison for the
# same reason - so both SL arms stay on the same single crossfit and differ in
# the pseudo-values alone.
pseudo_sl_t_standard <- function(
  X,
  pseudo,
  W,
  fold_indices,
  fold_list,
  sl_library = DEFAULT_SL_LIBRARY
) {
  n_obs <- nrow(X)
  cvps <- is.matrix(pseudo)

  tau_result <- future_map(
    seq_along(fold_list),
    function(i) {
      fold <- fold_list[i]
      in_train <- fold_indices != fold
      in_fold <- !in_train

      X_train <- X[in_train, , drop = FALSE]
      W_train <- W[in_train]
      pseudo_train <- if (cvps) pseudo[in_train, fold] else pseudo[in_train]
      X_test <- X[in_fold, , drop = FALSE]

      X0 <- as.data.frame(X_train[W_train == 0, , drop = FALSE])
      X1 <- as.data.frame(X_train[W_train == 1, , drop = FALSE])
      Y0 <- pseudo_train[W_train == 0]
      Y1 <- pseudo_train[W_train == 1]

      sl0 <- SuperLearner(
        Y = Y0,
        X = X0,
        SL.library = pretest_superlearner(Y0, X0, sl_library, gaussian()),
        cvControl = list(V = 5)
      )
      sl1 <- SuperLearner(
        Y = Y1,
        X = X1,
        SL.library = pretest_superlearner(Y1, X1, sl_library, gaussian()),
        cvControl = list(V = 5)
      )

      pred0 <- predict(sl0, newdata = as.data.frame(X_test))$pred
      pred1 <- predict(sl1, newdata = as.data.frame(X_test))$pred

      list(fold = fold, tau = as.numeric(pred1 - pred0))
    },
    .options = furrr_options(seed = TRUE)
  )

  tau <- rep(NA_real_, n_obs)
  for (result in tau_result) {
    tau[fold_indices == result$fold] <- result$tau
  }
  return(tau)
}

# SuperLearner T-learner: split pseudo-obs (Algorithm 2, Cwiling et al. 2025)
# 3-way split: training (V-2 folds), KM set (fold k+1 mod n_folds), validation (fold k)
# Training pseudo-obs: leave-one-out on training set; validation pseudo-obs: computed via KM set
pseudo_sl_t_split <- function(
  X,
  Y,
  D,
  W,
  horizon,
  fold_indices,
  fold_list,
  sl_library = DEFAULT_SL_LIBRARY
) {
  n_obs <- nrow(X)
  n_folds <- length(fold_list)
  D_int <- as.integer(D)
  Dc <- as.integer(D %in% c(1, 2))

  tau_result <- future_map(
    seq_along(fold_list),
    function(i) {
      fold <- fold_list[i]
      km_fold <- fold_list[(i %% n_folds) + 1] # next fold, wraps around

      in_val <- fold_indices == fold
      in_km <- fold_indices == km_fold
      in_train <- !in_val & !in_km

      X_train <- X[in_train, , drop = FALSE]
      W_train <- W[in_train]
      X_test <- X[in_val, , drop = FALSE]

      # Training pseudo-obs: standard leave-one-out on training set
      ps_RMTL_train <- pseudoyl(Y[in_train], D_int[in_train], horizon)
      ps_RMSTc_train <- pseudomean(Y[in_train], Dc[in_train], horizon)
      ps_RMTL1_train <- ps_RMTL_train$pseudo$cause1
      ps_RMTL2_train <- ps_RMTL_train$pseudo$cause2

      make_t_cate <- function(pseudo_train) {
        sl0 <- SuperLearner(
          Y = pseudo_train[W_train == 0],
          X = as.data.frame(X_train[W_train == 0, , drop = FALSE]),
          SL.library = sl_library,
          cvControl = list(V = 5)
        )
        sl1 <- SuperLearner(
          Y = pseudo_train[W_train == 1],
          X = as.data.frame(X_train[W_train == 1, , drop = FALSE]),
          SL.library = sl_library,
          cvControl = list(V = 5)
        )
        p0 <- predict(sl0, newdata = as.data.frame(X_test))$pred
        p1 <- predict(sl1, newdata = as.data.frame(X_test))$pred
        as.numeric(p1 - p0)
      }

      list(
        fold = fold,
        tau_RMTL1 = make_t_cate(ps_RMTL1_train),
        tau_RMTL2 = make_t_cate(ps_RMTL2_train),
        tau_RMSTc = make_t_cate(ps_RMSTc_train)
      )
    },
    .options = furrr_options(seed = TRUE)
  )

  tau_RMTL1 <- tau_RMTL2 <- tau_RMSTc <- rep(NA_real_, n_obs)
  for (result in tau_result) {
    idx <- fold_indices == result$fold
    tau_RMTL1[idx] <- result$tau_RMTL1
    tau_RMTL2[idx] <- result$tau_RMTL2
    tau_RMSTc[idx] <- result$tau_RMSTc
  }
  list(RMTL1 = tau_RMTL1, RMTL2 = tau_RMTL2, RMSTc = tau_RMSTc)
}

# SuperLearner DR-learner nuisances, single leave-one-fold-out ("scf_scf").
#
# Mirrors nuisance_pseudo_rf_scf but replaces regression_forest with
# SuperLearner, on the shape of R/cate_models.R::nuisance_sl. `pseudo` is the
# whole-sample vector (sl_dr_whole) or the n x V crossfit matrix (sl_dr_cvps);
# `pseudo_whole` remains the correction-term outcome in both, for the reason
# given on nuisance_pseudo_rf_scf.
nuisance_pseudo_sl <- function(
  X,
  pseudo,
  W,
  pseudo_whole,
  fold_indices,
  fold_list,
  sl_library = DEFAULT_SL_LIBRARY
) {
  cvps <- is.matrix(pseudo)

  cross_fits <- future_map(
    seq_along(fold_list),
    function(i) {
      fold <- fold_list[i]
      in_train <- fold_indices != fold
      in_test <- !in_train

      X_train <- as.data.frame(X[in_train, , drop = FALSE])
      X_test <- as.data.frame(X[in_test, , drop = FALSE])
      W_train <- W[in_train]
      W_test <- W[in_test]
      pseudo_train <- if (cvps) pseudo[in_train, fold] else pseudo[in_train]

      X_W_train <- cbind(W = W_train, X_train)

      sl_ps <- SuperLearner(
        Y = pseudo_train,
        X = X_W_train,
        SL.library = pretest_superlearner(
          pseudo_train,
          X_W_train,
          sl_library,
          gaussian()
        ),
        cvControl = list(V = 5)
      )
      sl_cf <- SuperLearner(
        Y = pseudo_train,
        X = X_train,
        SL.library = pretest_superlearner(
          pseudo_train,
          X_train,
          sl_library,
          gaussian()
        ),
        cvControl = list(V = 5)
      )
      sl_W <- SuperLearner(
        Y = W_train,
        X = X_train,
        SL.library = pretest_superlearner(
          W_train,
          X_train,
          sl_library,
          binomial()
        ),
        cvControl = list(V = 5)
      )

      pseudo0.hat <- as.numeric(
        predict(sl_ps, newdata = cbind(W = 0, X_test))$pred
      )
      pseudo1.hat <- as.numeric(
        predict(sl_ps, newdata = cbind(W = 1, X_test))$pred
      )
      pseudo.hat.cf <- as.numeric(predict(sl_cf, newdata = X_test)$pred)
      W.hat <- as.numeric(predict(sl_W, newdata = X_test)$pred)

      # failsafes if SuperLearner returns all-zero predictions, as R/cate_models.R
      if (all(pseudo0.hat == 0) && all(pseudo1.hat == 0)) {
        warning("SuperLearner failed for pseudo.hat. Using mean(pseudo).")
        pseudo0.hat <- rep(
          mean(pseudo_train[W_train == 0], na.rm = TRUE),
          sum(in_test)
        )
        pseudo1.hat <- rep(
          mean(pseudo_train[W_train == 1], na.rm = TRUE),
          sum(in_test)
        )
      }
      if (all(W.hat == 0)) {
        warning("SuperLearner failed for W.hat. Using mean(W).")
        W.hat <- rep(mean(W_train, na.rm = TRUE), sum(in_test))
      }

      W.hat <- trim_ps(W.hat)

      pseudo.hat <- W_test * pseudo1.hat + (1 - W_test) * pseudo0.hat
      po <- dr_pseudo(
        pseudo_whole[in_test],
        W_test,
        pseudo1.hat,
        pseudo0.hat,
        W.hat
      )

      list(
        fold = fold,
        po = po,
        pseudo.hat = pseudo.hat,
        pseudo0.hat = pseudo0.hat,
        pseudo.hat.cf = pseudo.hat.cf,
        W.hat = W.hat
      )
    },
    .options = furrr_options(seed = TRUE)
  )

  scatter_folds(
    cross_fits,
    fold_indices,
    c("po", "pseudo.hat", "pseudo0.hat", "pseudo.hat.cf", "W.hat")
  )
}

# Stage 2 is R/cate_models.R::stage_2_sl, whose is.vector(po) branch handles the
# vector `po` these nuisances now return and routes it through
# pretest_superlearner. The local matrix-branch copy this study used to carry is
# gone with the double crossfitting that produced the matrix.
#
# X is coerced to a data.frame here because stage_2_sl indexes it straight into
# SuperLearner(), which warns ("X is not a data frame") and can silently drop
# candidate learners on a matrix. R/cate_models.R does the same coercion at its
# own SuperLearner call site (cate_methods, `X <- as.data.frame(X)`); this study
# keeps a matrix X for the grf arms, so the conversion belongs here.
pseudo_dr_sl <- function(
  X,
  po,
  fold_indices,
  fold_list,
  sl_library = DEFAULT_SL_LIBRARY
) {
  stage_2_sl(as.data.frame(X), po, fold_indices, fold_list, sl_library)
}

# SuperLearner DR-learner nuisances: split pseudo-obs
# 3-way split per fold: training (V-2 folds), KM set (fold k+1), validation (fold k)
# Validation pseudo-obs are split pseudo-obs (independent of training), per Algorithm 2
nuisance_pseudo_sl_split <- function(
  X,
  Y,
  D,
  W,
  horizon,
  fold_indices,
  fold_list,
  sl_library = DEFAULT_SL_LIBRARY
) {
  n_obs <- nrow(X)
  n_folds <- length(fold_list)
  D_int <- as.integer(D)
  Dc <- as.integer(D %in% c(1, 2))

  cross_fits <- future_map(
    seq_along(fold_list),
    function(i) {
      fold <- fold_list[i]
      km_fold <- fold_list[(i %% n_folds) + 1]

      in_val <- fold_indices == fold
      in_km <- fold_indices == km_fold
      in_train <- !in_val & !in_km

      X_train <- as.data.frame(X[in_train, , drop = FALSE])
      X_val <- as.data.frame(X[in_val, , drop = FALSE])
      W_train <- W[in_train]
      W_val <- W[in_val]

      # Standard pseudo-obs on training set
      ps_RMTL_train <- pseudoyl(Y[in_train], D_int[in_train], horizon)
      ps_RMSTc_train <- pseudomean(Y[in_train], Dc[in_train], horizon)

      # Split pseudo-obs for validation fold using KM set
      split_RMTL <- compute_split_pseudoyl(
        Y[in_km],
        D_int[in_km],
        Y[in_val],
        D_int[in_val],
        horizon
      )
      split_RMSTc <- compute_split_pseudomean(
        Y[in_km],
        Dc[in_km],
        Y[in_val],
        Dc[in_val],
        horizon
      )

      make_dr_po <- function(pseudo_train, pseudo_val_split) {
        sl_ps <- SuperLearner(
          Y = pseudo_train,
          X = cbind(W = W_train, X_train),
          SL.library = sl_library,
          cvControl = list(V = 5)
        )
        sl_W <- SuperLearner(
          Y = W_train,
          X = X_train,
          SL.library = sl_library,
          cvControl = list(V = 5)
        )

        pseudo0.hat <- as.numeric(
          predict(sl_ps, newdata = cbind(W = 0, X_val))$pred
        )
        pseudo1.hat <- as.numeric(
          predict(sl_ps, newdata = cbind(W = 1, X_val))$pred
        )
        W.hat <- as.numeric(predict(sl_W, newdata = X_val)$pred)

        cate <- pseudo1.hat - pseudo0.hat
        pseudo.hat <- W_val * pseudo1.hat + (1 - W_val) * pseudo0.hat
        po <- cate +
          ((pseudo_val_split - pseudo.hat) * (W_val - W.hat)) /
            (W.hat * (1 - W.hat))
        po
      }

      list(
        fold = fold,
        in_val = which(in_val),
        po_RMTL1 = make_dr_po(ps_RMTL_train$pseudo$cause1, split_RMTL$cause1),
        po_RMTL2 = make_dr_po(ps_RMTL_train$pseudo$cause2, split_RMTL$cause2),
        po_RMSTc = make_dr_po(ps_RMSTc_train, split_RMSTc)
      )
    },
    .options = furrr_options(seed = TRUE)
  )

  po_RMTL1 <- po_RMTL2 <- po_RMSTc <- rep(NA_real_, n_obs)
  for (result in cross_fits) {
    po_RMTL1[result$in_val] <- result$po_RMTL1
    po_RMTL2[result$in_val] <- result$po_RMTL2
    po_RMSTc[result$in_val] <- result$po_RMSTc
  }
  list(po_RMTL1 = po_RMTL1, po_RMTL2 = po_RMTL2, po_RMSTc = po_RMSTc)
}

# SuperLearner stage 2 for split DR-learner (po is an n-vector, not a matrix)
#
# Now functionally redundant with R/cate_models.R::stage_2_sl's is.vector(po)
# branch, which additionally pretests the library. Kept because the disabled
# sl_dr_split block in all_cate_surv_models() and surv_dr_split_na_diagnose.R
# both reference it by name, and the split-pseudo arms are out of scope for the
# crossfitting change. Fold it into stage_2_sl when that NA bug is fixed.
stage_2_sl_vec <- function(
  X,
  po_vec,
  fold_indices,
  fold_list,
  sl_library = DEFAULT_SL_LIBRARY
) {
  n_obs <- nrow(X)

  tau_results <- future_map(
    seq_along(fold_list),
    function(i) {
      fold <- fold_list[i]
      in_train <- fold_indices != fold
      in_fold <- !in_train

      sl <- SuperLearner(
        Y = po_vec[in_train],
        X = as.data.frame(X[in_train, , drop = FALSE]),
        SL.library = sl_library,
        cvControl = list(V = 5)
      )
      list(
        fold = fold,
        predictions = as.numeric(
          predict(sl, newdata = as.data.frame(X[in_fold, , drop = FALSE]))$pred
        )
      )
    },
    .options = furrr_options(seed = TRUE)
  )

  tau <- rep(NA_real_, n_obs)
  for (result in tau_results) {
    tau[fold_indices == result$fold] <- result$predictions
  }
  tau
}
