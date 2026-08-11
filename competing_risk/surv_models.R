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

  # fold creation for crossfitting
  fold_indices <- sort(seq(n_obs) %% n_folds) + 1
  fold_list <- unique(fold_indices)
  fold_pairs <- utils::combn(fold_list, 2, simplify = FALSE)

  # all results container
  results <- list()

  message("Causal Forest IPW approaches (RMST1, RMST2, RMSTc)...")
  results$ipw <- list()
  results$ipw$RMST1 <- cf_ipw(
    X,
    Y,
    D,
    W,
    horizon,
    fold_indices,
    fold_list,
    event = 1
  )
  results$ipw$RMST2 <- cf_ipw(
    X,
    Y,
    D,
    W,
    horizon,
    fold_indices,
    fold_list,
    event = 2
  )
  results$ipw$RMSTc <- cf_ipw(
    X,
    Y,
    D,
    W,
    horizon,
    fold_indices,
    fold_list,
    event = "composite"
  )

  message("Causal Survival Forest cs approaches (RMST1, RMST2, RMSTc)...")
  results$csf_cs <- list()
  results$csf_cs$RMST1 <- csf_cs(
    X,
    Y,
    D,
    W,
    horizon,
    fold_indices,
    fold_list,
    event = 1
  )
  results$csf_cs$RMST2 <- csf_cs(
    X,
    Y,
    D,
    W,
    horizon,
    fold_indices,
    fold_list,
    event = 2
  )
  results$csf_cs$RMSTc <- csf_cs(
    X,
    Y,
    D,
    W,
    horizon,
    fold_indices,
    fold_list,
    event = "composite"
  )

  message("Causal Survival Forest sh approaches (RMST1, RMST2)...")
  results$csf_sh <- list()
  results$csf_sh$RMST1 <- csf_sh(
    X,
    Y,
    D,
    W,
    horizon,
    fold_indices,
    fold_list,
    event = 1
  )
  results$csf_sh$RMST2 <- csf_sh(
    X,
    Y,
    D,
    W,
    horizon,
    fold_indices,
    fold_list,
    event = 2
  )

  message("Pseudo-value approaches...")
  pseudo_whole <- pseudo_all(Y, D, horizon)

  message("Pseudo-value Causal Forest (RMTL1, RMTL2, RMSTc)...")
  results$pseudo_cf <- list()
  pseudo_cv <- pseudo_crossfit(
    Y,
    D,
    horizon,
    fold_indices,
    fold_list,
    pseudo_whole
  )
  results$pseudo_cf$RMTL1 <- pseudo_cf(
    X,
    pseudo_cv$ps_RMTL1,
    W,
    fold_indices,
    fold_list
  )
  results$pseudo_cf$RMTL2 <- pseudo_cf(
    X,
    pseudo_cv$ps_RMTL2,
    W,
    fold_indices,
    fold_list
  )
  results$pseudo_cf$RMSTc <- pseudo_cf(
    X,
    pseudo_cv$ps_RMSTc,
    W,
    fold_indices,
    fold_list
  )

  message("Pseudo-value DR Learner...")
  results$pseudo_dr <- list()
  message("Jack-knife pseudovalues for the DR learner...") # I've written pseudo so much I'm convinced it's not a word anymore
  pseudo_dr_cv <- pseudo_double_crossfit(
    Y,
    D,
    horizon,
    fold_indices,
    fold_pairs,
    pseudo_whole
  )
  message("DR nuisance function estimation (RMTL1, RMTL2, RMSTc)...")
  nuisance_RMTL1 <- nuisance_pseudo_rf(
    X,
    pseudo_dr_cv$ps_RMTL1,
    W,
    pseudo_whole$ps_RMTL1,
    fold_indices,
    fold_pairs
  )
  nuisance_RMTL2 <- nuisance_pseudo_rf(
    X,
    pseudo_dr_cv$ps_RMTL2,
    W,
    pseudo_whole$ps_RMTL2,
    fold_indices,
    fold_pairs
  )
  nuisance_RMSTc <- nuisance_pseudo_rf(
    X,
    pseudo_dr_cv$ps_RMSTc,
    W,
    pseudo_whole$ps_RMSTc,
    fold_indices,
    fold_pairs
  )

  message(" DR Second stage regression (RMTL1, RMTL2, RMSTc)...")
  results$pseudo_dr$RMTL1 <- pseudo_dr_rf(
    X,
    nuisance_RMTL1$po_matrix,
    fold_indices,
    fold_list
  )
  results$pseudo_dr$RMTL2 <- pseudo_dr_rf(
    X,
    nuisance_RMTL2$po_matrix,
    fold_indices,
    fold_list
  )
  results$pseudo_dr$RMSTc <- pseudo_dr_rf(
    X,
    nuisance_RMSTc$po_matrix,
    fold_indices,
    fold_list
  )

  message("SuperLearner T-learner (standard pseudo-obs)...")
  results$sl_t_std <- list()
  results$sl_t_std$RMTL1 <- pseudo_sl_t_standard(
    X,
    pseudo_cv$ps_RMTL1,
    W,
    fold_indices,
    fold_list,
    sl_library
  )
  results$sl_t_std$RMTL2 <- pseudo_sl_t_standard(
    X,
    pseudo_cv$ps_RMTL2,
    W,
    fold_indices,
    fold_list,
    sl_library
  )
  results$sl_t_std$RMSTc <- pseudo_sl_t_standard(
    X,
    pseudo_cv$ps_RMSTc,
    W,
    fold_indices,
    fold_list,
    sl_library
  )

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

  message("SuperLearner DR-learner (standard pseudo-obs)...")
  results$sl_dr_std <- list()
  nuisance_sl_RMTL1 <- nuisance_pseudo_sl(
    X,
    pseudo_dr_cv$ps_RMTL1,
    W,
    pseudo_whole$ps_RMTL1,
    fold_indices,
    fold_pairs,
    sl_library
  )
  nuisance_sl_RMTL2 <- nuisance_pseudo_sl(
    X,
    pseudo_dr_cv$ps_RMTL2,
    W,
    pseudo_whole$ps_RMTL2,
    fold_indices,
    fold_pairs,
    sl_library
  )
  nuisance_sl_RMSTc <- nuisance_pseudo_sl(
    X,
    pseudo_dr_cv$ps_RMSTc,
    W,
    pseudo_whole$ps_RMSTc,
    fold_indices,
    fold_pairs,
    sl_library
  )
  results$sl_dr_std$RMTL1 <- pseudo_dr_sl(
    X,
    nuisance_sl_RMTL1$po_matrix,
    fold_indices,
    fold_list,
    sl_library
  )
  results$sl_dr_std$RMTL2 <- pseudo_dr_sl(
    X,
    nuisance_sl_RMTL2$po_matrix,
    fold_indices,
    fold_list,
    sl_library
  )
  results$sl_dr_std$RMSTc <- pseudo_dr_sl(
    X,
    nuisance_sl_RMSTc$po_matrix,
    fold_indices,
    fold_list,
    sl_library
  )

  # SuperLearner DR-learner (split pseudo-obs) - DISABLED for now, see
  # README.md "Known issues": compute_split_pseudoyl()/compute_split_pseudomean()
  # (below) hit an unpatched NaN in pseudo::pseudoyl()'s internal ci.omit()
  # whenever a validation Y exceeds the max Y in that fold's KM set, and
  # unlike pseudo_crossfit()/pseudo_double_crossfit() there is no NA
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
  results$pseudos$dr_cv <- pseudo_dr_cv
  results$pseudos$cf_cv <- pseudo_cv
  results$pseudos$whole <- pseudo_whole

  results$nuisances <- list()
  results$nuisances$RMTL1 <- nuisance_RMTL1
  results$nuisances$RMTL2 <- nuisance_RMTL2
  results$nuisances$RMSTc <- nuisance_RMSTc

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
cf_ipw <- function(X, Y, D, W, horizon, fold_indices, fold_list, event = 1) {
  n_obs <- nrow(X)

  RMST_result <- future_map(
    seq_along(fold_list),
    function(i) {
      fold <- fold_list[i]
      in_train <- fold_indices != fold
      in_fold <- !in_train

      # get ipw weights
      weights_0 <- get_ipw(
        X[in_train, ],
        Y[in_train],
        D[in_train],
        W[in_train],
        horizon,
        0
      ) # when there is no censoring these are all 1's

      if (event == "composite") {
        total_observed <- weights_0$observed
        sample_weights <- weights_0$ipw[total_observed]
      } else {
        # Identify event of interest
        all_events <- unique(D[D != 0])
        competing <- setdiff(all_events, event)

        weights_competing <- get_ipw(
          X[in_train, ],
          Y[in_train],
          D[in_train],
          W[in_train],
          horizon,
          competing
        ) # weights to account for the competing event

        total_observed <- weights_0$observed & weights_competing$observed
        sample_weights <- weights_0$ipw * weights_competing$ipw
        sample_weights <- sample_weights[total_observed]
      }

      # keep non censored/CE obs in training folds
      train_idx <- which(in_train)
      include <- train_idx[total_observed]

      # causal forest (not survival because we've sorted out weights - see grf tutorial)
      # Y truncated at horizon so the outcome is min(T1*, horizon), matching the RMST estimand
      forest <- causal_forest(
        X[include, ],
        pmin(Y[include], horizon),
        W[include],
        sample.weights = sample_weights
      )

      pred <- predict(forest, newdata = X[in_fold, ])

      list(fold = fold, tau = pred$predictions)
    },
    .options = furrr_options(seed = TRUE)
  )

  # put predictions into the right order
  tau_RMST <- rep(NA_real_, n_obs)
  for (result in RMST_result) {
    in_fold <- fold_indices == result$fold
    tau_RMST[in_fold] <- result$tau
  }

  return(tau_RMST)
}
# CSF - treating CEs as censoring events
csf_cs <- function(X, Y, D, W, horizon, fold_indices, fold_list, event = 1) {
  n_obs <- nrow(X)

  if (event == "composite") {
    # Modify event to include 1 and 2
    D_event <- as.numeric(D %in% c(1, 2))
  } else {
    # Modify event def to treat competing event as censoring
    D_event <- as.numeric(D == event)
  }

  RMST_result <- future_map(
    seq_along(fold_list),
    function(i) {
      fold <- fold_list[i]
      in_train <- fold_indices != fold
      in_fold <- !in_train

      forest <- causal_survival_forest(
        X[in_train, ],
        Y[in_train],
        W[in_train],
        D_event[in_train],
        target = "RMST",
        horizon = horizon
      )

      pred <- predict(forest, newdata = X[in_fold, ])

      list(fold = fold, tau = pred$predictions)
    },
    .options = furrr_options(seed = TRUE)
  )

  # put predictions back
  tau_RMST <- rep(NA, n_obs)
  for (result in RMST_result) {
    in_fold <- fold_indices == result$fold
    tau_RMST[in_fold] <- result$tau
  }

  return(tau_RMST)
}
# CSF - keep CEs in the risk set
csf_sh <- function(X, Y, D, W, horizon, fold_indices, fold_list, event = 1) {
  n_obs <- nrow(X)

  # move competing events after horizon (keep them in the risk set)
  D_sh <- as.numeric(D == event)
  Y_sh <- ifelse(!(D %in% c(event, 0)), horizon + 1, Y)

  RMST_result <- future_map(
    seq_along(fold_list),
    function(i) {
      fold <- fold_list[i]
      in_train <- fold_indices != fold
      in_fold <- !in_train

      include <- in_train
      # if there is censoring for event horizon, account for this
      cens <- any(D[in_train] == 0 & Y[in_train] < horizon)
      if (cens) {
        weights_0 <- get_ipw(
          X[in_train, ],
          Y[in_train],
          D[in_train],
          W[in_train],
          horizon,
          0
        ) # when there is no censoring these are all 1's

        total_observed <- weights_0$observed
        sample_weights <- weights_0$ipw[total_observed]
        train_idx <- which(in_train)
        include <- train_idx[total_observed]
      }

      forest <- causal_survival_forest(
        X[include, ],
        Y_sh[include],
        W[include],
        D_sh[include],
        target = "RMST",
        horizon = horizon,
        sample.weights = if (cens) sample_weights else NULL
      )

      pred <- predict(forest, newdata = X[in_fold, ])

      list(fold = fold, tau = pred$predictions)
    },
    .options = furrr_options(seed = TRUE)
  )

  # put predictions back
  tau_RMST <- rep(NA, n_obs)
  for (result in RMST_result) {
    in_fold <- fold_indices == result$fold
    tau_RMST[in_fold] <- result$tau
  }

  return(tau_RMST)
}
# pseudovalue aproaches
pseudo_all <- function(Y, D, horizon) {
  # reformatting for different estimands
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
# crossfit pseudovalues to use in a causal forest?
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
  Y_sh <- ifelse(D == 2, horizon + 1, Y)
  D_sh <- as.integer(D == 1)

  pseudos <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold <- !in_train

    ps_RMTL <- pseudoyl(Y[in_train], D_int[in_train], horizon)
    ps_RMSTc <- pseudomean(Y[in_train], Dc[in_train], horizon)
    ps_sh_RMST <- pseudomean(Y_sh[in_train], D_sh[in_train], horizon)

    # check for any failed computations of the values, replace with value from pseudo_whole
    # this is not a perfect solution - introduces some data leakage ?
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

  # Empty matrices for pseudos
  ps_RMTL1_mat <- matrix(NA_real_, nrow = n_obs, ncol = n_folds)
  ps_RMTL2_mat <- matrix(NA_real_, nrow = n_obs, ncol = n_folds)
  ps_RMSTc_mat <- matrix(NA_real_, nrow = n_obs, ncol = n_folds)
  ps_sh_RMST_mat <- matrix(NA_real_, nrow = n_obs, ncol = n_folds)

  # Fill matrices
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
# double crossfitting for pseudo outcomes to use in DR learner
pseudo_double_crossfit <- function(
  Y,
  D,
  horizon,
  fold_indices,
  fold_pairs,
  pseudo_whole
) {
  n_obs <- length(Y)
  n_pairs <- length(fold_pairs)
  # reformatting for different estimands
  D_int <- as.integer(D)
  Dc <- as.integer(D %in% c(1, 2))
  Y_sh <- ifelse(D == 2, horizon + 1, Y)
  D_sh <- as.integer(D == 1)

  pseudos <- future_map(seq_along(fold_pairs), function(i) {
    fold_pair <- fold_pairs[[i]]
    in_train <- !(fold_indices %in% fold_pair)
    in_test <- !in_train

    ps_RMTL <- pseudoyl(Y[in_train], D_int[in_train], horizon)
    ps_RMSTc <- pseudomean(Y[in_train], Dc[in_train], horizon)
    ps_sh_RMST <- pseudomean(Y_sh[in_train], D_sh[in_train], horizon)

    # check for any failed computations of the values, replace with value from pseudo_whole
    # this is not a perfect solution - introduces some data leakage ?
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

  # make matrices
  fold_list <- unique(fold_indices)
  ps_RMTL1_mat <- matrix(NA_real_, nrow = n_obs, ncol = length(fold_pairs))
  ps_RMTL2_mat <- matrix(NA_real_, nrow = n_obs, ncol = length(fold_pairs))
  ps_RMSTc_mat <- matrix(NA_real_, nrow = n_obs, ncol = length(fold_pairs))
  ps_sh_RMST_mat <- matrix(NA_real_, nrow = n_obs, ncol = length(fold_pairs))

  for (i in seq_along(fold_pairs)) {
    fold_pair <- fold_pairs[[i]]
    in_train <- !(fold_indices %in% fold_pair)

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
# causal forest using pseudo values
pseudo_cf <- function(X, pseudo, W, fold_indices, fold_list) {
  n_obs <- nrow(X)

  tau_result <- future_map(
    seq_along(fold_list),
    function(i) {
      fold <- fold_list[i]
      in_train <- fold_indices != fold
      in_fold <- !in_train

      forest <- causal_forest(
        X[in_train, ],
        pseudo[in_train, fold],
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
nuisance_pseudo_rf <- function(
  X,
  pseudo,
  W,
  pseudo_whole,
  fold_indices,
  fold_pairs
) {
  cross_fits <- future_map(
    seq_along(fold_pairs),
    function(i) {
      fold_pair <- fold_pairs[[i]]
      in_train <- !(fold_indices %in% fold_pair)
      in_test <- !in_train

      ps.hat.model <- regression_forest(
        cbind(W[in_train], X[in_train, ]),
        pseudo[in_train, i]
      )
      ps.hat.cf.model <- regression_forest(X[in_train, ], pseudo[in_train, i])
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
      W.hat <- predict(W.hat.model, newdata = X_test)$predictions

      # DR learner pseudo outcome (not the same as the pseudo values we already have)
      W_test <- W[in_test]
      pseudo.hat <- W_test * pseudo1.hat + (1 - W_test) * pseudo0.hat

      cate <- pseudo1.hat - pseudo0.hat

      # use pseudo values fit on the whole data set as the Y's in the pseudo outcome equation?
      pseudo_test <- pseudo_whole[in_test]
      po <- cate +
        ((pseudo_test - pseudo.hat) * (W_test - W.hat)) / (W.hat * (1 - W.hat))

      list(
        po = po,
        pseudo.hat = pseudo.hat,
        pseudo0.hat = pseudo0.hat,
        pseudo.hat.cf = pseudo.hat.cf,
        W.hat = W.hat,
        fold_pair = fold_pair
      )
    },
    .options = furrr_options(seed = TRUE)
  )

  # matrices of nuisances
  fold_list <- unique(fold_indices)
  po_matrix <- collate_predictions(
    fold_list,
    fold_pairs,
    fold_indices,
    cross_fits,
    "po"
  )
  pseudo.hat_matrix <- collate_predictions(
    fold_list,
    fold_pairs,
    fold_indices,
    cross_fits,
    "pseudo.hat"
  )
  pseudo.hat.cf_matrix <- collate_predictions(
    fold_list,
    fold_pairs,
    fold_indices,
    cross_fits,
    "pseudo.hat.cf"
  )
  pseudo0.hat_matrix <- collate_predictions(
    fold_list,
    fold_pairs,
    fold_indices,
    cross_fits,
    "pseudo0.hat"
  )
  W.hat_matrix <- collate_predictions(
    fold_list,
    fold_pairs,
    fold_indices,
    cross_fits,
    "W.hat"
  )

  list(
    po_matrix = po_matrix,
    po = rowMeans(po_matrix, na.rm = T),
    pseudo.hat_matrix = pseudo.hat_matrix,
    pseudo.hat = rowMeans(pseudo.hat_matrix, na.rm = T),
    pseudo.hat.cf_matrix = pseudo.hat.cf_matrix,
    pseudo.hat.cf = rowMeans(pseudo.hat.cf_matrix, na.rm = T),
    pseudo0.hat_matrix = pseudo0.hat_matrix,
    pseudo0.hat = rowMeans(pseudo0.hat_matrix, na.rm = T),
    W.hat_matrix = W.hat_matrix,
    W.hat = rowMeans(W.hat_matrix, na.rm = T)
  )
}

stage_2_rf <- function(X, po, fold_indices, fold_list) {
  n_obs <- nrow(X)

  tau_results <- future_map(
    seq_along(fold_list),
    function(i) {
      fold <- fold_list[i]
      in_train <- fold_indices != fold
      in_fold <- !in_train

      forest <- regression_forest(X[in_train, ], po[in_train, fold])

      tau_pred <- predict(forest, newdata = X[in_fold, ])$predictions

      list(fold = fold, predictions = tau_pred)
    },
    .options = furrr_options(seed = TRUE)
  )

  # Reconstruct tau vector
  tau <- rep(NA, n_obs)
  for (result in tau_results) {
    fold <- result$fold
    in_fold <- fold_indices == fold
    tau[in_fold] <- result$predictions
  }
  return(tau)
}

pseudo_dr_rf <- function(X, po, fold_indices, fold_list) {
  tau <- stage_2_rf(X, po, fold_indices, fold_list)

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

# SuperLearner T-learner: standard cross-fit pseudo-obs
# pseudo_cv: n x n_folds matrix from pseudo_crossfit (leave-fold-out)
pseudo_sl_t_standard <- function(
  X,
  pseudo_cv,
  W,
  fold_indices,
  fold_list,
  sl_library = DEFAULT_SL_LIBRARY
) {
  n_obs <- nrow(X)

  tau_result <- future_map(
    seq_along(fold_list),
    function(i) {
      fold <- fold_list[i]
      in_train <- fold_indices != fold
      in_fold <- !in_train

      X_train <- X[in_train, , drop = FALSE]
      W_train <- W[in_train]
      pseudo_train <- pseudo_cv[in_train, i]
      X_test <- X[in_fold, , drop = FALSE]

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

# SuperLearner DR-learner nuisances: standard cross-fit pseudo-obs
# Mirrors nuisance_pseudo_rf but replaces regression_forest with SuperLearner
nuisance_pseudo_sl <- function(
  X,
  pseudo,
  W,
  pseudo_whole,
  fold_indices,
  fold_pairs,
  sl_library = DEFAULT_SL_LIBRARY
) {
  cross_fits <- future_map(
    seq_along(fold_pairs),
    function(i) {
      fold_pair <- fold_pairs[[i]]
      in_train <- !(fold_indices %in% fold_pair)
      in_test <- !in_train

      X_train <- as.data.frame(X[in_train, , drop = FALSE])
      X_test <- as.data.frame(X[in_test, , drop = FALSE])
      W_train <- W[in_train]
      W_test <- W[in_test]

      sl_ps <- SuperLearner(
        Y = pseudo[in_train, i],
        X = cbind(W = W_train, X_train),
        SL.library = sl_library,
        cvControl = list(V = 5)
      )
      sl_cf <- SuperLearner(
        Y = pseudo[in_train, i],
        X = X_train,
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
        predict(sl_ps, newdata = cbind(W = 0, X_test))$pred
      )
      pseudo1.hat <- as.numeric(
        predict(sl_ps, newdata = cbind(W = 1, X_test))$pred
      )
      pseudo.hat.cf <- as.numeric(predict(sl_cf, newdata = X_test)$pred)
      W.hat <- as.numeric(predict(sl_W, newdata = X_test)$pred)

      pseudo.hat <- W_test * pseudo1.hat + (1 - W_test) * pseudo0.hat
      cate <- pseudo1.hat - pseudo0.hat
      pseudo_test <- pseudo_whole[in_test]
      po <- cate +
        ((pseudo_test - pseudo.hat) * (W_test - W.hat)) / (W.hat * (1 - W.hat))

      list(
        po = po,
        pseudo.hat = pseudo.hat,
        pseudo0.hat = pseudo0.hat,
        pseudo.hat.cf = pseudo.hat.cf,
        W.hat = W.hat,
        fold_pair = fold_pair
      )
    },
    .options = furrr_options(seed = TRUE)
  )

  fold_list <- unique(fold_indices)
  po_matrix <- collate_predictions(
    fold_list,
    fold_pairs,
    fold_indices,
    cross_fits,
    "po"
  )
  pseudo.hat_matrix <- collate_predictions(
    fold_list,
    fold_pairs,
    fold_indices,
    cross_fits,
    "pseudo.hat"
  )
  pseudo0.hat_matrix <- collate_predictions(
    fold_list,
    fold_pairs,
    fold_indices,
    cross_fits,
    "pseudo0.hat"
  )
  W.hat_matrix <- collate_predictions(
    fold_list,
    fold_pairs,
    fold_indices,
    cross_fits,
    "W.hat"
  )

  list(
    po_matrix = po_matrix,
    po = rowMeans(po_matrix, na.rm = TRUE),
    pseudo.hat_matrix = pseudo.hat_matrix,
    pseudo.hat = rowMeans(pseudo.hat_matrix, na.rm = TRUE),
    pseudo0.hat_matrix = pseudo0.hat_matrix,
    pseudo0.hat = rowMeans(pseudo0.hat_matrix, na.rm = TRUE),
    W.hat_matrix = W.hat_matrix,
    W.hat = rowMeans(W.hat_matrix, na.rm = TRUE)
  )
}

# SuperLearner stage 2 for DR-learner (po is an n x n_folds matrix)
stage_2_sl <- function(
  X,
  po,
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
  return(tau)
}

pseudo_dr_sl <- function(
  X,
  po,
  fold_indices,
  fold_list,
  sl_library = DEFAULT_SL_LIBRARY
) {
  stage_2_sl(X, po, fold_indices, fold_list, sl_library)
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
