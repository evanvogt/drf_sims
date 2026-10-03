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
library(randomForestSRC)

# sl_libraries, sl_fit_predict, pretest_superlearner. Also arrives via
# R/cate_models.R, but sourced here too so the diagnostic scripts that load
# this file alone still find it.
source(here::here("R", "sl_library.R"))

# Per-nuisance libraries. surv_analysis.R passes sl_libraries(n); this default
# is the full, n > 100 set. Until the library change it was the single
# vector c("SL.glm", "SL.glmnet", "SL.ranger", "SL.gam") for every nuisance.
DEFAULT_SL_LIBRARY <- sl_libraries(Inf)

# num.threads is grf's thread count, forwarded to every grf fit and predict()
# below, as R/cate_models.R::cate_methods does. NULL (default) is grf's own
# default (all visible cores), so callers that omit it - the diagnostic scripts,
# R/regression_check.R - behave as before. surv_analysis.R passes grf_threads.
# The SuperLearner arms are single-threaded already (SL.ranger's num.threads = 1).
# rfsrc has no thread argument: its fold workers set rf.cores = 1 (see
# nuisance_rsf_scf), and its main-process OOB fits follow OMP_NUM_THREADS, which
# surv_analysis.R sets to grf_threads.
all_cate_surv_models <- function(
  data,
  n_folds = 10,
  horizon = 30,
  sl_library = DEFAULT_SL_LIBRARY,
  num.threads = NULL
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
  results$ipw$RMST1 <- cf_ipw(X, Y, D, W, horizon, event = 1, num.threads)
  results$ipw$RMST2 <- cf_ipw(X, Y, D, W, horizon, event = 2, num.threads)
  results$ipw$RMSTc <- cf_ipw(X, Y, D, W, horizon, event = "composite", num.threads)

  message("Causal Survival Forest cs approaches (RMST1, RMST2, RMSTc)...")
  results$csf_cs <- list()
  results$csf_cs$RMST1 <- csf_cs(X, Y, D, W, horizon, event = 1, num.threads)
  results$csf_cs$RMST2 <- csf_cs(X, Y, D, W, horizon, event = 2, num.threads)
  results$csf_cs$RMSTc <- csf_cs(X, Y, D, W, horizon, event = "composite", num.threads)

  message("Causal Survival Forest sh approaches (RMST1, RMST2)...")
  results$csf_sh <- list()
  results$csf_sh$RMST1 <- csf_sh(X, Y, D, W, horizon, event = 1, num.threads)
  results$csf_sh$RMST2 <- csf_sh(X, Y, D, W, horizon, event = 2, num.threads)

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
    pseudo_cf_whole_oob(X, ps_whole(e), W, num.threads)
  })
  results$pseudo_cf_whole_scf <- by_estimand(function(e) {
    pseudo_cf_scf(X, ps_whole(e), W, fold_indices, fold_list, num.threads)
  })
  results$pseudo_cf_cvps_scf <- by_estimand(function(e) {
    pseudo_cf_scf(X, ps_cv(e), W, fold_indices, fold_list, num.threads)
  })

  message("Pseudo-value DR Learner nuisances (whole_oob, whole_scf, cvps_scf)...")
  nuis_dr <- list(
    whole_oob = by_estimand(function(e) {
      nuisance_pseudo_rf_oob(X, ps_whole(e), W, num.threads)
    }),
    whole_scf = by_estimand(function(e) {
      nuisance_pseudo_rf_scf(
        X,
        ps_whole(e),
        W,
        ps_whole(e),
        fold_indices,
        fold_list,
        num.threads
      )
    }),
    cvps_scf = by_estimand(function(e) {
      nuisance_pseudo_rf_scf(
        X,
        ps_cv(e),
        W,
        ps_whole(e),
        fold_indices,
        fold_list,
        num.threads
      )
    })
  )

  message("DR second stage regression (RMTL1, RMTL2, RMSTc)...")
  results$pseudo_dr_whole_oob <- by_estimand(function(e) {
    stage2_whole_rf(X, nuis_dr$whole_oob[[e]]$po, num.threads = num.threads)$tau
  })
  results$pseudo_dr_whole_scf <- by_estimand(function(e) {
    stage_2_rf_scf(X, nuis_dr$whole_scf[[e]]$po, fold_indices, fold_list, num.threads)
  })
  results$pseudo_dr_cvps_scf <- by_estimand(function(e) {
    stage_2_rf_scf(X, nuis_dr$cvps_scf[[e]]$po, fold_indices, fold_list, num.threads)
  })

  message("SuperLearner T-learner (whole and cvps pseudo-obs)...")
  results$sl_t_whole <- by_estimand(function(e) {
    pseudo_sl_t_standard(X, ps_whole(e), W, fold_indices, fold_list, sl_library)
  })
  results$sl_t_cvps <- by_estimand(function(e) {
    pseudo_sl_t_standard(X, ps_cv(e), W, fold_indices, fold_list, sl_library)
  })

  message("SuperLearner T-learner (split pseudo-obs)...")
  # pseudo_whole is passed for its NA fallback only - see pseudo_sl_t_split. It
  # does NOT make this a "whole pseudo-values" arm: the training pseudo-values
  # are still computed on the training folds, and results$sl_t_split$n_na_fallback
  # counts how many rows had to borrow from the whole-sample vector.
  results$sl_t_split <- pseudo_sl_t_split(
    X,
    Y,
    D,
    W,
    horizon,
    fold_indices,
    fold_list,
    sl_library,
    pseudo_whole = pseudo_whole
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

  # Last of the fitted arms on purpose: rfsrc() draws its seed from R's RNG, so
  # placing it here leaves every arm above on the stream it had before.
  message("Random survival forest DR-learner (oob, scf)...")
  nuis_rsf <- list(
    oob = nuisance_rsf_oob(X, Y, D, W, horizon, pseudo_whole, num.threads),
    scf = nuisance_rsf_scf(
      X,
      Y,
      D,
      W,
      horizon,
      pseudo_whole,
      fold_indices,
      fold_list,
      num.threads
    )
  )
  results$rsf_dr_oob <- by_estimand(function(e) {
    stage2_whole_rf(X, nuis_rsf$oob[[e]]$po, num.threads = num.threads)$tau
  })
  results$rsf_dr_scf <- by_estimand(function(e) {
    stage_2_rf_scf(X, nuis_rsf$scf[[e]]$po, fold_indices, fold_list, num.threads)
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
  results$nuisances <- list(rf = nuis_dr, sl = nuis_sl, rsf = nuis_rsf)

  results$fold_indices <- fold_indices

  return(results)
}

# Helper Functions
# get IPW weights for specific event (censoring or the competing event)
get_ipw <- function(X, Y, D, W, horizon, censor, num.threads = NULL) {
  # select censoring event
  D_event <- as.numeric(D == censor)

  # fit survival forest for time to censoring event
  sf_censor <- survival_forest(cbind(X, W), Y, D_event, num.threads = num.threads)

  censor_prob <- predict(
    sf_censor,
    failure.times = pmin(Y, horizon),
    prediction.times = "time",
    num.threads = num.threads
  )$predictions

  # get weights for non censored (clipped)
  observed <- (D != censor)
  epsilon <- 1e-3
  ipw <- 1 / pmax(censor_prob, epsilon)

  list(observed = observed, ipw = ipw)
}
# CF using IPW weights for censoring and competing events - see GRF tutorial
#
# Whole-sample, grf-internal crossfitting ("cf_default"), matching
# R/cate_models.R::run_causal_forest. The forest is fit only on `include` -
# observations not censored and free of the competing event - so OOB predictions
# exist for those rows only; the excluded rows never entered the forest at all,
# so a plain newdata prediction for them is honest. Same pattern as
# crossfitting/cf_models.R::nuisance_oob_rf.
cf_ipw <- function(X, Y, D, W, horizon, event = 1, num.threads = NULL) {
  n_obs <- nrow(X)

  # get ipw weights. get_ipw()'s predict() call passes no newdata, so the
  # censoring probabilities are already grf OOB predictions.
  weights_0 <- get_ipw(X, Y, D, W, horizon, 0, num.threads) # all 1's when there is no censoring

  if (event == "composite") {
    total_observed <- weights_0$observed
    sample_weights <- weights_0$ipw[total_observed]
  } else {
    # Identify the competing event (the cause not targeted)
    all_events <- unique(D[D != 0])
    competing <- setdiff(all_events, event)

    weights_competing <- get_ipw(X, Y, D, W, horizon, competing, num.threads)

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
    sample.weights = sample_weights,
    num.threads = num.threads
  )

  tau_RMST <- rep(NA_real_, n_obs)
  tau_RMST[include] <- predict(forest, num.threads = num.threads)$predictions # OOB for the rows it was fit on
  if (length(include) < n_obs) {
    tau_RMST[-include] <- predict(
      forest,
      newdata = X[-include, , drop = FALSE],
      num.threads = num.threads
    )$predictions
  }

  return(tau_RMST)
}
# CSF - treating competing events as censoring events
# Whole sample; causal_survival_forest cross-fits its own nuisances internally
# and predict() with no newdata returns OOB predictions.
csf_cs <- function(X, Y, D, W, horizon, event = 1, num.threads = NULL) {
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
    horizon = horizon,
    num.threads = num.threads
  )

  return(predict(forest, num.threads = num.threads)$predictions)
}
# CSF - keep competing events in the risk set
# Whole sample, as csf_cs. When censoring is present the forest is fit on the
# uncensored subset only, so the excluded rows get a newdata prediction - see
# cf_ipw above for why that is honest.
csf_sh <- function(X, Y, D, W, horizon, event = 1, num.threads = NULL) {
  n_obs <- nrow(X)

  # move competing events after horizon (keep them in the risk set)
  D_sh <- as.numeric(D == event)
  Y_sh <- ifelse(!(D %in% c(event, 0)), horizon + 1, Y)

  include <- seq_len(n_obs)
  sample_weights <- NULL
  # if there is censoring for event horizon, account for this
  cens <- any(D == 0 & Y < horizon)
  if (cens) {
    weights_0 <- get_ipw(X, Y, D, W, horizon, 0, num.threads)
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
    sample.weights = sample_weights,
    num.threads = num.threads
  )

  tau_RMST <- rep(NA_real_, n_obs)
  tau_RMST[include] <- predict(forest, num.threads = num.threads)$predictions
  if (length(include) < n_obs) {
    tau_RMST[-include] <- predict(
      forest,
      newdata = X[-include, , drop = FALSE],
      num.threads = num.threads
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
pseudo_cf_whole_oob <- function(X, pseudo, W, num.threads = NULL) {
  forest <- causal_forest(X, pseudo, W, num.threads = num.threads)
  predict(forest, num.threads = num.threads)$predictions
}

# causal forest using pseudo values - single leave-one-fold-out ("scf").
#
# `pseudo` is either the whole-sample vector (the whole_scf control arm) or the
# n x V crossfit matrix from pseudo_crossfit (the cvps_scf arm). Branching on
# is.matrix() here is what lets one function serve both arms, so that the two
# differ in the pseudo-values alone and in nothing else.
pseudo_cf_scf <- function(X, pseudo, W, fold_indices, fold_list, num.threads = NULL) {
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
        W[in_train],
        num.threads = num.threads
      )

      pred <- predict(forest, newdata = X[in_fold, ], num.threads = num.threads)

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
# Whole-sample OOB, T-learner ("oob_oob") - the pseudo-value counterpart of
# R/cate_models.R::nuisance_rf, and the production-parity arm. No sample
# splitting: t_learner_rf (one forest per arm; own-arm predictions OOB) supplies
# both counterfactuals, and two more whole-sample forests supply W.hat and
# pseudo.hat.cf, all taken out-of-bag. Until the move to per-arm outcome models
# this was the S-learner "oob_oob_s" (one forest on cbind(W, X), read through
# grf's X.orig).
#
# `pseudo` is the whole-sample pseudo-value vector, and it is also the outcome in
# the DR correction term - there is no separate `pseudo_whole` argument here
# because with no split the two coincide.
nuisance_pseudo_rf_oob <- function(X, pseudo, W, num.threads = NULL) {
  mu <- t_learner_rf(X, pseudo, W, num.threads = num.threads)
  pseudo0.hat <- mu$Y0.hat
  pseudo1.hat <- mu$Y1.hat

  W.hat <- trim_ps(predict(regression_forest(X, W, num.threads = num.threads))$predictions)
  pseudo.hat.cf <- predict(regression_forest(X, pseudo, num.threads = num.threads))$predictions

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
  fold_list,
  num.threads = NULL
) {
  cvps <- is.matrix(pseudo)

  cross_fits <- future_map(
    seq_along(fold_list),
    function(i) {
      fold <- fold_list[i]
      in_train <- fold_indices != fold
      in_test <- !in_train

      pseudo_train <- if (cvps) pseudo[in_train, fold] else pseudo[in_train]
      X_train <- X[in_train, ]
      W_train <- W[in_train]

      # one outcome forest per arm (T-learner), control arm first
      arm_forest <- function(arm) {
        regression_forest(X_train[W_train == arm, , drop = FALSE],
                          pseudo_train[W_train == arm],
                          num.threads = num.threads)
      }
      ps0.model <- arm_forest(0)
      ps1.model <- arm_forest(1)
      ps.hat.cf.model <- regression_forest(X_train, pseudo_train, num.threads = num.threads)
      W.hat.model <- regression_forest(X_train, W_train, num.threads = num.threads)

      X_test <- X[in_test, ]

      pred_test <- function(model) {
        predict(model, newdata = X_test, num.threads = num.threads)$predictions
      }
      pseudo0.hat <- pred_test(ps0.model)
      pseudo1.hat <- pred_test(ps1.model)
      pseudo.hat.cf <- pred_test(ps.hat.cf.model)
      W.hat <- trim_ps(pred_test(W.hat.model))

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
stage_2_rf_scf <- function(X, po, fold_indices, fold_list, num.threads = NULL) {
  n_obs <- nrow(X)

  tau_results <- future_map(
    seq_along(fold_list),
    function(i) {
      fold <- fold_list[i]
      in_train <- fold_indices != fold
      in_fold <- !in_train

      forest <- regression_forest(X[in_train, ], po[in_train], num.threads = num.threads)

      tau_pred <- predict(forest, newdata = X[in_fold, ], num.threads = num.threads)$predictions

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

# Random survival forest DR-learner (randomForestSRC)
#
# The outcome model is a competing-risks random survival forest fit to the
# observed (Y, D), not a regression on pseudo-values: one forest per arm
# (T-learner), default composite Gray splitting (splitrule = "logrankCR"), so a
# single fit gives both CIFs and, from them, all three estimands. Pseudo-values
# enter only the DR correction term, as whole-sample values, so the whole/cvps
# factor the grf arms vary has nothing to act on here - rsf_dr_oob and
# rsf_dr_scf differ in fitting alone. Valid because censoring in this DGM is
# independent of X, so E[pseudo | X, W] is the conditional RMTL the forest
# estimates. W.hat and stage 2 are the grf ones the pseudo_dr arms use, so the
# outcome model is the only thing that differs from pseudo_dr_whole_*.
RSF_ESTIMANDS <- c("RMTL1", "RMTL2", "RMSTc")
RSF_NUISANCES <- c("po", "pseudo.hat", "pseudo0.hat", "pseudo.hat.cf", "W.hat")

#' Competing-risks forest on (Y, D), censored at the horizon
#'
#' Censoring at the horizon leaves F_j(t) on [0, horizon] unchanged and keeps
#' Gray's split statistic to the estimand's window. ntime = 0 keeps every event
#' time in time.interest (the default, 150, thins it to a grid), which
#' rsf_rmtl's tail term relies on. rfsrc() draws its seed from R's runif() when
#' none is given, so it follows setup_rng_stream() and furrr's seeds as the grf
#' fits do. The event codes present are kept as attr(, "causes") for rsf_rmtl's
#' single-cause fallback, since a predict() object does not carry them.
rsf_fit <- function(X, Y, D, horizon) {
  status <- as.integer(ifelse(Y > horizon, 0L, D))
  df <- data.frame(time = pmin(Y, horizon), status = status, X)
  fit <- rfsrc(Surv(time, status) ~ ., data = df, ntime = 0)
  attr(fit, "causes") <- sort(unique(status[status > 0]))
  fit
}

#' RMTL1, RMTL2 and RMSTc to the horizon from an rfsrc grow or predict object
#'
#' For family "surv-CR", predicted / predicted.oob is the package's "expected
#' number of life years lost due to cause j": sum_{q<T} CIF_j(t_q) (t_{q+1} -
#' t_q) over time.interest (src/survival.c, getMortality), the step-function CIF
#' integrated from the first event time to the LAST one, t_T. That upper limit
#' cannot be set - ntime only snaps to observed event times - so the piece from
#' t_T to the horizon, (horizon - t_T) * CIF_j(t_T), is added here. Nothing is
#' missing below t_1, where the CIF is 0. The package's own get.rmst() is not
#' used: it is unexported, "surv" family only, and takes the in-bag survival
#' whenever the OOB one exists. CIF columns are in sorted event-code order.
#'
#' @param o rfsrc grow object (oob = TRUE) or predict() object (oob = FALSE)
#' @param causes event codes in the training rows - attr(fit, "causes")
#' @return n x 3 matrix, columns RSF_ESTIMANDS
rsf_rmtl <- function(o, horizon, oob, causes) {
  times <- o$time.interest
  n_t <- length(times)
  tail <- horizon - times[n_t]

  if (o$family == "surv-CR") {
    years_lost <- if (oob) o$predicted.oob else o$predicted
    cif <- if (oob) o$cif.oob else o$cif
    rmtl <- years_lost + tail * matrix(cif[, n_t, ], ncol = dim(cif)[3])
  } else {
    # Only one cause in the training rows: rfsrc falls back to a plain survival
    # forest, whose `predicted` is ensemble mortality, not years lost. 1 - S is
    # then the present cause's CIF; the absent cause's AJ estimate is 0.
    warning("rsf_rmtl: only cause ", paste(causes, collapse = ", "),
            " in the training rows; the other cause's RMTL is set to 0.")
    surv <- if (oob) o$survival.oob else o$survival
    surv <- matrix(surv, ncol = n_t)
    rmtl <- matrix(0, nrow(surv), 2)
    rmtl[, causes] <- (1 - surv) %*% c(diff(times), tail)
  }

  cbind(RMTL1 = rmtl[, 1], RMTL2 = rmtl[, 2],
        RMSTc = horizon - rmtl[, 1] - rmtl[, 2])
}

#' Predict RMTLs at new rows from a fitted rsf_fit() forest
rsf_predict_rmtl <- function(fit, X_new, horizon) {
  pred <- predict(fit, newdata = as.data.frame(X_new))
  rsf_rmtl(pred, horizon, oob = FALSE, causes = attr(fit, "causes"))
}

#' The DR nuisance list for each estimand, from n x 3 outcome-model matrices
#'
#' Same fields as nuisance_pseudo_rf_oob / _scf return, so
#' surv_nuisance_extract.R reads these arms unchanged. `pseudo` is
#' pseudo_all()'s list (ps_RMTL1, ps_RMTL2, ps_RMSTc), already subset to the
#' rows in hand.
rsf_dr_nuisances <- function(pseudo, W, mu0, mu1, mu_cf, W.hat) {
  setNames(lapply(RSF_ESTIMANDS, function(e) {
    list(
      po = dr_pseudo(pseudo[[paste0("ps_", e)]], W, mu1[, e], mu0[, e], W.hat),
      pseudo.hat = W * mu1[, e] + (1 - W) * mu0[, e],
      pseudo0.hat = mu0[, e],
      pseudo.hat.cf = mu_cf[, e],
      W.hat = W.hat
    )
  }), RSF_ESTIMANDS)
}

#' Whole-sample OOB RSF nuisances ("oob")
#'
#' The rfsrc counterpart of t_learner_rf (R/cate_models.R): one forest per arm,
#' each unit's own-arm prediction from predicted.oob / cif.oob and its other-arm
#' prediction from a forest that never saw it. pseudo.hat.cf comes from a pooled
#' forest on X alone (OOB), as nuisance_pseudo_rf_oob's does from grf.
nuisance_rsf_oob <- function(X, Y, D, W, horizon, pseudo_whole, num.threads = NULL) {
  n_obs <- nrow(X)

  arm_fit <- function(arm) {
    in_arm <- W == arm
    fit <- rsf_fit(X[in_arm, , drop = FALSE], Y[in_arm], D[in_arm], horizon)
    pred <- matrix(NA_real_, n_obs, length(RSF_ESTIMANDS),
                   dimnames = list(NULL, RSF_ESTIMANDS))
    pred[in_arm, ] <- rsf_rmtl(fit, horizon, oob = TRUE,
                               causes = attr(fit, "causes"))
    pred[!in_arm, ] <- rsf_predict_rmtl(fit, X[!in_arm, , drop = FALSE], horizon)
    pred
  }
  # control arm first, as t_learner_rf
  mu0 <- arm_fit(0)
  mu1 <- arm_fit(1)

  cf_fit <- rsf_fit(X, Y, D, horizon)
  mu_cf <- rsf_rmtl(cf_fit, horizon, oob = TRUE, causes = attr(cf_fit, "causes"))

  W.hat <- trim_ps(predict(regression_forest(X, W, num.threads = num.threads))$predictions)

  rsf_dr_nuisances(pseudo_whole, W, mu0, mu1, mu_cf, W.hat)
}

#' Single leave-one-fold-out RSF nuisances ("scf")
#'
#' The crossfit twin of nuisance_rsf_oob, on the shape of nuisance_pseudo_rf_scf:
#' every forest is fit on the training folds and predicts the held-out fold.
#' Whole-sample pseudo-values stay in the correction term, as in every DR arm.
nuisance_rsf_scf <- function(
  X,
  Y,
  D,
  W,
  horizon,
  pseudo_whole,
  fold_indices,
  fold_list,
  num.threads = NULL
) {
  cross_fits <- future_map(
    seq_along(fold_list),
    function(i) {
      # rfsrc's OpenMP would otherwise take every core in each fold worker
      op <- options(rf.cores = 1)
      on.exit(options(op), add = TRUE)

      fold <- fold_list[i]
      in_train <- fold_indices != fold
      in_test <- !in_train

      X_train <- X[in_train, , drop = FALSE]
      X_test <- X[in_test, , drop = FALSE]
      W_train <- W[in_train]

      fit_predict <- function(rows) {
        fit <- rsf_fit(X_train[rows, , drop = FALSE], Y[in_train][rows],
                       D[in_train][rows], horizon)
        rsf_predict_rmtl(fit, X_test, horizon)
      }
      # control arm first, as nuisance_rsf_oob
      mu0 <- fit_predict(W_train == 0)
      mu1 <- fit_predict(W_train == 1)
      mu_cf <- fit_predict(rep(TRUE, sum(in_train)))

      W.hat.model <- regression_forest(X_train, W_train, num.threads = num.threads)
      W.hat <- trim_ps(predict(W.hat.model, newdata = X_test,
                               num.threads = num.threads)$predictions)

      list(
        fold = fold,
        nuis = rsf_dr_nuisances(
          lapply(pseudo_whole, `[`, in_test),
          W[in_test],
          mu0,
          mu1,
          mu_cf,
          W.hat
        )
      )
    },
    .options = furrr_options(seed = TRUE)
  )

  setNames(lapply(RSF_ESTIMANDS, function(e) {
    per_fold <- lapply(cross_fits, function(r) c(list(fold = r$fold), r$nuis[[e]]))
    scatter_folds(per_fold, fold_indices, RSF_NUISANCES)
  }), RSF_ESTIMANDS)
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

# One pseudoyl() call per validation unit returns both causes; this used to
# make two identical calls, one per cause.
compute_split_pseudoyl <- function(Y_km, D_km, Y_val, D_val, horizon) {
  n_km <- length(Y_km)
  ps <- vapply(
    seq_along(Y_val),
    function(j) {
      p <- pseudoyl(
        c(Y_km, Y_val[j]),
        as.integer(c(D_km, D_val[j])),
        horizon
      )$pseudo
      c(p$cause1[n_km + 1], p$cause2[n_km + 1])
    },
    numeric(2)
  )
  list(cause1 = ps[1, ], cause2 = ps[2, ])
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
  # one outcome regression per arm - the per-arm outcome library
  arm_lib <- as_sl_libs(sl_library)$Y

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
      newX <- list(test = as.data.frame(X_test))

      pred0 <- sl_fit_predict(
        Y0,
        X0,
        newX,
        pretest_superlearner(Y0, X0, arm_lib, gaussian()),
        cvControl = list(V = 5)
      )$pred$test
      pred1 <- sl_fit_predict(
        Y1,
        X1,
        newX,
        pretest_superlearner(Y1, X1, arm_lib, gaussian()),
        cvControl = list(V = 5)
      )$pred$test

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
#
# This function used to take down the whole run, and was the ONLY arm that did:
# 225 of the 1400 array indices never produced a res_sim_*.RDS, and
# surv_failed_diagnose.R traced every reproducible one of them to here, by two
# routes. It was the only SuperLearner call site in this study that both passed
# `sl_library` unvalidated AND computed its pseudo-values with an unguarded
# pseudoyl(), so it was missing both of the guards the rest of the file has:
#
#   * NA pseudo-values (195 of the 225). pseudoyl()/pseudomean() return NA for
#     whoever holds the maximum observed time in the sample being decomposed -
#     the pseudo::ci.omit bug written up in README "Known issues". The NA went
#     straight into SuperLearner(), which refuses it:
#       "missing data is currently not supported. Check Y, X, and newX"
#     Now guarded with the same is.na() -> pseudo_whole substitution
#     pseudo_crossfit() has used all along, and counted the same way - see
#     n_na_fallback below.
#
#   * a degenerate library (25 of the 225). Where the treatment effect on the
#     event 1 is strong (bW_1 = -0.7, scenarios 1/3/4/6/7, under the
#     pre-2026-10-01 parameters), nearly
#     every treated subject has the cause-1 event before the horizon, so the
#     treated arm's RMTL2 pseudo-values collapse onto a handful of distinct
#     values - 3 across 195 rows on array index 67. SL.glmnet then dies inside
#     SuperLearner's own CV ("y is constant; gaussian glmnet fails at
#     standardization step"), leaving NULL in $fitLibrary, and predict()'s
#     default onlySL = FALSE predicted from that NULL fit:
#       "no applicable method for 'predict' applied to an object of class NULL"
#     Now pretested, exactly as pseudo_sl_t_standard() already was.
#
# @param pseudo_whole whole-sample pseudo-values from pseudo_all(), used only to
#   fill NAs. Passed in rather than recomputed: all_cate_surv_models() already
#   has it. Optional so the function still stands alone for diagnostics; without
#   it the NA rows are dropped from the fold's fit instead.
pseudo_sl_t_split <- function(
  X,
  Y,
  D,
  W,
  horizon,
  fold_indices,
  fold_list,
  sl_library = DEFAULT_SL_LIBRARY,
  pseudo_whole = NULL
) {
  n_obs <- nrow(X)
  n_folds <- length(fold_list)
  D_int <- as.integer(D)
  Dc <- as.integer(D %in% c(1, 2))
  # one outcome regression per arm - the per-arm outcome library
  arm_lib <- as_sl_libs(sl_library)$Y

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
      X_test <- as.data.frame(X[in_val, , drop = FALSE])

      # Training pseudo-obs: standard leave-one-out on training set
      ps_RMTL_train <- pseudoyl(Y[in_train], D_int[in_train], horizon)
      ps_RMSTc_train <- pseudomean(Y[in_train], Dc[in_train], horizon)
      ps_RMTL1_train <- ps_RMTL_train$pseudo$cause1
      ps_RMTL2_train <- ps_RMTL_train$pseudo$cause2

      # Count before substituting, so the leakage the substitution introduces
      # stays measurable - the same contract pseudo_crossfit() keeps. A
      # whole-sample pseudo-value has seen the validation and KM folds, so a
      # non-zero count qualifies this arm's independence claim.
      n_na <- c(
        RMTL1 = sum(is.na(ps_RMTL1_train)),
        RMTL2 = sum(is.na(ps_RMTL2_train)),
        RMSTc = sum(is.na(ps_RMSTc_train))
      )

      fill_na <- function(ps, whole_name) {
        if (!anyNA(ps)) return(ps)
        if (is.null(pseudo_whole)) return(ps) # fit_arm's `keep` drops them
        ifelse(is.na(ps), pseudo_whole[[whole_name]][in_train], ps)
      }
      ps_RMTL1_train <- fill_na(ps_RMTL1_train, "ps_RMTL1")
      ps_RMTL2_train <- fill_na(ps_RMTL2_train, "ps_RMTL2")
      ps_RMSTc_train <- fill_na(ps_RMSTc_train, "ps_RMSTc")

      # One treatment arm's fit. `keep` is the backstop for the two cases
      # fill_na cannot cover: no pseudo_whole passed at all, and the rare one
      # where the whole-sample vector is itself NA on that row (pseudoyl fails
      # on the GLOBAL max-time observation too). Either way no NA reaches
      # SuperLearner, which is what aborted 195 of the 225 runs.
      fit_arm <- function(pseudo_train, w) {
        in_arm <- W_train == w
        y <- pseudo_train[in_arm]
        x <- as.data.frame(X_train[in_arm, , drop = FALSE])
        keep <- !is.na(y)

        lib <- pretest_superlearner(y[keep], x[keep, , drop = FALSE],
                                    arm_lib, gaussian())
        # predicted at newX during the fit (sl_fit_predict), from the weighted
        # learners only, so a candidate that survives pretest's 2-fold CV and
        # then fails in the live 5-fold one cannot reintroduce the NULL-fit
        # crash above. (This was predict(onlySL = TRUE) before the fit moved to
        # sl_fit_predict.)
        pred <- sl_fit_predict(
          y[keep],
          x[keep, , drop = FALSE],
          list(test = X_test),
          lib,
          cvControl = list(V = 5)
        )$pred$test

        # Same failsafe as R/cate_models.R::nuisance_sl: when every learner ends
        # up with zero weight SuperLearner returns all-zero predictions, which
        # are not an estimate of anything. Testing for EXACT zeros is safe even
        # for RMTL2, where the treated arm's pseudo-values genuinely sit near
        # zero - a real fit does not return floating-point 0 on every row.
        # anyNA is folded in so a stray NA falls back rather than propagating
        # into tau, which is the failure this whole function was fixed for.
        if (anyNA(pred) || isTRUE(all(pred == 0))) {
          warning("pseudo_sl_t_split: SuperLearner gave no usable predictions ",
                  "for W = ", w, " on fold ", fold, ". Using mean(pseudo).")
          pred <- rep(mean(y[keep], na.rm = TRUE), nrow(X_test))
        }
        pred
      }

      # Control arm first, as before the fix - these consume the RNG stream, so
      # keeping the order avoids a gratuitous change to the numbers on top of
      # the one pretest_superlearner already makes.
      make_t_cate <- function(pseudo_train) {
        p0 <- fit_arm(pseudo_train, 0)
        p1 <- fit_arm(pseudo_train, 1)
        p1 - p0
      }

      list(
        fold = fold,
        n_na = n_na,
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
  list(
    RMTL1 = tau_RMTL1,
    RMTL2 = tau_RMTL2,
    RMSTc = tau_RMSTc,
    n_na_fallback = colSums(do.call(rbind, lapply(tau_result, `[[`, "n_na")))
  )
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
  libs <- as_sl_libs(sl_library)

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

      # one outcome model per arm (T-learner), control arm first
      arm_pred <- function(arm) {
        y <- pseudo_train[W_train == arm]
        x <- X_train[W_train == arm, , drop = FALSE]
        sl_fit_predict(
          y,
          x,
          list(test = X_test),
          pretest_superlearner(y, x, libs$Y, gaussian()),
          cvControl = list(V = 5)
        )$pred$test
      }
      pseudo0.hat <- arm_pred(0)
      pseudo1.hat <- arm_pred(1)

      # the marginal outcome model, on every training row
      sl_cf <- sl_fit_predict(
        pseudo_train,
        X_train,
        list(test = X_test),
        pretest_superlearner(
          pseudo_train,
          X_train,
          libs$Y,
          gaussian()
        ),
        cvControl = list(V = 5)
      )
      # binomial, as it is pretested and as R/cate_models.R::nuisance_sl fits
      # it. Until the library change this fit took SuperLearner's default
      # gaussian family (NNLS), after a binomial pretest.
      sl_W <- sl_fit_predict(
        W_train,
        X_train,
        list(test = X_test),
        pretest_superlearner(W_train, X_train, libs$W, binomial()),
        family = binomial(),
        cvControl = list(V = 5)
      )

      pseudo.hat.cf <- sl_cf$pred$test
      W.hat <- sl_W$pred$test

      # failsafes if SuperLearner returns all-zero predictions, as
      # R/cate_models.R - per arm, since each arm now has its own model
      if (all(pseudo0.hat == 0)) {
        warning("SuperLearner failed for pseudo.hat in arm W = 0. Using its mean.")
        pseudo0.hat <- rep(mean(pseudo_train[W_train == 0], na.rm = TRUE), sum(in_test))
      }
      if (all(pseudo1.hat == 0)) {
        warning("SuperLearner failed for pseudo.hat in arm W = 1. Using its mean.")
        pseudo1.hat <- rep(mean(pseudo_train[W_train == 1], na.rm = TRUE), sum(in_test))
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
  tau <- stage_2_sl(as.data.frame(X), po, fold_indices, fold_list, sl_library)
  # results here are one tau vector per estimand; the dropped-learner table
  # stage_2_sl attaches is printed by the pretest as it runs
  attr(tau, "sl_dropped") <- NULL
  tau
}

# SuperLearner DR-learner nuisances: split pseudo-obs
# DISABLED arm (see all_cate_surv_models).
# 3-way split per fold: training (V-2 folds), KM set (fold k+1), validation (fold k)
# Validation pseudo-obs are split pseudo-obs (independent of training), per Algorithm 2
#
# Nuisances as nuisance_pseudo_sl: one outcome model per arm (T-learner,
# pretested) and a pretested binomial propensity, trimmed by trim_ps. Until
# 2026-10-01 the outcome model was an S-learner on cbind(W, X), the propensity
# was untrimmed, and it was refit once per estimand on identical data; it is
# now fit once per fold. The all-zero failsafes use isTRUE() so that an NA
# training pseudo-value (pseudoyl's max-time bug) still flows through to NA in
# po rather than erroring here - surv_dr_split_na_diagnose.R relies on that.
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
  libs <- as_sl_libs(sl_library)

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

      # propensity: does not depend on the estimand, so fit once per fold
      W.hat <- sl_fit_predict(
        W_train,
        X_train,
        list(val = X_val),
        pretest_superlearner(W_train, X_train, libs$W, binomial()),
        family = binomial(),
        cvControl = list(V = 5)
      )$pred$val
      if (isTRUE(all(W.hat == 0))) {
        warning("SuperLearner failed for W.hat. Using mean(W).")
        W.hat <- rep(mean(W_train, na.rm = TRUE), sum(in_val))
      }
      W.hat <- trim_ps(W.hat)

      make_dr_po <- function(pseudo_train, pseudo_val_split) {
        # one outcome model per arm (T-learner), control arm first
        arm_pred <- function(arm) {
          y <- pseudo_train[W_train == arm]
          x <- X_train[W_train == arm, , drop = FALSE]
          pred <- sl_fit_predict(
            y,
            x,
            list(val = X_val),
            pretest_superlearner(y, x, libs$Y, gaussian()),
            cvControl = list(V = 5)
          )$pred$val
          if (isTRUE(all(pred == 0))) {
            warning("SuperLearner failed for pseudo.hat in arm W = ", arm,
                    ". Using its mean.")
            pred <- rep(mean(y, na.rm = TRUE), sum(in_val))
          }
          pred
        }
        pseudo0.hat <- arm_pred(0)
        pseudo1.hat <- arm_pred(1)

        dr_pseudo(pseudo_val_split, W_val, pseudo1.hat, pseudo0.hat, W.hat)
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

      sl <- sl_fit_predict(
        po_vec[in_train],
        as.data.frame(X[in_train, , drop = FALSE]),
        list(fold = as.data.frame(X[in_fold, , drop = FALSE])),
        as_sl_libs(sl_library)$tau,
        cvControl = list(V = 5)
      )
      list(fold = fold, predictions = sl$pred$fold)
    },
    .options = furrr_options(seed = TRUE)
  )

  tau <- rep(NA_real_, n_obs)
  for (result in tau_results) {
    tau[fold_indices == result$fold] <- result$predictions
  }
  tau
}
