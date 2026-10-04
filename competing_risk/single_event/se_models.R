##########
# Title: running CATE models - single-event survival
##########
# A subset of competing_risk/surv_models.R's arms, on one event. Expects
# R/cate_models.R and competing_risk/surv_models.R to be sourced first, in that
# order (as se_analysis.R does), and reuses from them unchanged:
#   csf_cs, pseudo_cf_whole_oob, nuisance_pseudo_rf_oob, nuisance_pseudo_sl,
#   pseudo_dr_sl, stage_2_rf_scf (surv_models.R);
#   stage2_whole_rf, dr_pseudo (R/cate_models.R); scatter_folds, trim_ps (R/utils.R)
#
# What is new here:
#   - pseudo_rmst(): KM pseudo-values of the RMST, the single-event counterpart
#     of pseudo_all();
#   - the RSF DR-learner on a survival-family forest. The parent's rsf_fit /
#     rsf_rmtl handle one cause only through their single-cause warning branch,
#     which would fire on every fit here;
#   - t_from_dr(): the T-learners, read off the DR-learners' per-arm outcome
#     models at no extra fitting cost.
#
# The estimand is the RMST to the horizon. Every arm returns one tau vector,
# stored as results$<framework>$RMST.

#' Fit every CATE arm on one single-event dataset
#'
#' @param data data.frame as generate_se_data() emits it: Y, D (0/1), W, then X
#' @param sl_library SuperLearner libraries, as sl_libraries(n) builds them
#' @param num.threads grf's thread count; NULL is grf's default
all_cate_se_models <- function(
  data,
  n_folds = 10,
  horizon = 28,
  sl_library = DEFAULT_SL_LIBRARY,
  num.threads = NULL
) {
  X <- as.matrix(data[, !(names(data) %in% c("Y", "D", "W"))])
  Y <- data$Y
  D <- data$D
  W <- data$W
  n_obs <- nrow(X)

  # folds for the single-crossfit arms (SuperLearner, rsf_*_scf)
  fold_indices <- sort(seq(n_obs) %% n_folds) + 1
  fold_list <- unique(fold_indices)

  results <- list()

  message("Causal survival forest...")
  results$csf <- list(RMST = csf_cs(X, Y, D, W, horizon, event = 1, num.threads))

  message("Pseudo-values...")
  ps <- pseudo_rmst(Y, D, horizon)

  message("Pseudo-value causal forest...")
  results$pseudo_cf_whole_oob <- list(RMST = pseudo_cf_whole_oob(X, ps, W, num.threads))

  message("Pseudo-value RF DR- and T-learner...")
  nuis_rf <- nuisance_pseudo_rf_oob(X, ps, W, num.threads)
  results$pseudo_dr_whole_oob <- list(
    RMST = stage2_whole_rf(X, nuis_rf$po, num.threads = num.threads)$tau
  )
  results$pseudo_t_whole_oob <- list(RMST = t_from_dr(nuis_rf, ps, W))

  message("Pseudo-value SuperLearner DR- and T-learner...")
  nuis_sl <- nuisance_pseudo_sl(X, ps, W, ps, fold_indices, fold_list, sl_library)
  results$sl_dr_whole <- list(
    RMST = pseudo_dr_sl(X, nuis_sl$po, fold_indices, fold_list, sl_library)
  )
  results$sl_t_whole <- list(RMST = t_from_dr(nuis_sl, ps, W))

  # Last of the fitted arms on purpose: rfsrc() draws its seed from R's RNG, so
  # placing it here leaves every arm above on the stream it had before.
  message("Random survival forest DR- and T-learner (oob, scf)...")
  nuis_rsf <- list(
    oob = nuisance_rsf_se_oob(X, Y, D, W, horizon, ps, num.threads),
    scf = nuisance_rsf_se_scf(X, Y, D, W, horizon, ps, fold_indices, fold_list,
                              num.threads)
  )
  results$rsf_dr_oob <- list(
    RMST = stage2_whole_rf(X, nuis_rsf$oob$po, num.threads = num.threads)$tau
  )
  results$rsf_t_oob <- list(RMST = t_from_dr(nuis_rsf$oob, ps, W))
  results$rsf_dr_scf <- list(
    RMST = stage_2_rf_scf(X, nuis_rsf$scf$po, fold_indices, fold_list, num.threads)
  )
  results$rsf_t_scf <- list(RMST = t_from_dr(nuis_rsf$scf, ps, W))

  results$pseudos <- list(whole = ps)
  results$nuisances <- list(rf = nuis_rf, sl = nuis_sl, rsf = nuis_rsf)
  results$fold_indices <- fold_indices

  results
}

#' Jackknife pseudo-values of the RMST from one Kaplan-Meier fit on all n rows
#'
#' No covariates are needed: censoring is independent of X and W. An NA
#' pseudo-value is an error, not a fallback - competing_risk's route C (README
#' "Known issues") came from an unguarded NA here. It needs nobody to be
#' observed past the horizon, which this DGM makes vanishingly unlikely (about
#' 28% of controls are event-free at 28), so a run that does hit it should fail
#' and show up in check_failed() rather than be patched over.
pseudo_rmst <- function(Y, D, horizon) {
  ps <- pseudomean(Y, as.integer(D), horizon)
  if (anyNA(ps)) {
    stop("pseudo_rmst: ", sum(is.na(ps)), " NA pseudo-value(s); max(Y) = ",
         signif(max(Y), 4), ", horizon = ", horizon)
  }
  as.numeric(ps)
}

#' T-learner CATE from a DR-learner's nuisances
#'
#' dr_pseudo() returns po = mu1 - mu0 + (theta - mu_W)(W - e) / (e (1 - e)),
#' and every nuisance list saves pseudo.hat = mu_W and the trimmed W.hat it used,
#' so mu1 - mu0 is po minus that correction - exact algebra, no refit. The
#' T-learner is then the same per-arm outcome models the DR-learner used, with
#' the same honesty (OOB own-arm predictions, other-arm ones from a model that
#' never saw the unit, or single crossfit).
#'
#' @param nuis a nuisance list with po, pseudo.hat and W.hat
#' @param theta the pseudo-values in the DR correction term
t_from_dr <- function(nuis, theta, W) {
  nuis$po - (theta - nuis$pseudo.hat) * (W - nuis$W.hat) /
    (nuis$W.hat * (1 - nuis$W.hat))
}

# ---- random survival forest DR-learner (randomForestSRC, survival family) ----
#
# As competing_risk's rsf_dr_* arms, on one event: one survival forest per arm
# (T-learner), default log-rank splitting, mu = the forest's survival curve
# integrated to the horizon. Pseudo-values enter only the DR correction term.
# W.hat and stage 2 are the grf ones the pseudo_dr arm uses, so the outcome
# model is the only thing that differs from pseudo_dr_whole_oob.

#' Survival forest on (Y, D), censored at the horizon
#'
#' Censoring at the horizon leaves S(t) on [0, horizon] unchanged. ntime = 0
#' keeps every event time in time.interest, which rsf_rmst's tail term relies on.
rsf_se_fit <- function(X, Y, D, horizon) {
  status <- as.integer(ifelse(Y > horizon, 0L, D))
  df <- data.frame(time = pmin(Y, horizon), status = status, X)
  rfsrc(Surv(time, status) ~ ., data = df, ntime = 0)
}

#' RMST to the horizon from an rfsrc grow (oob = TRUE) or predict object
#'
#' S is a step function on time.interest, equal to 1 before the first event
#' time, so integral_0^h S = h - integral (1 - S) = h - sum_q (1 - S(t_q)) *
#' (t_{q+1} - t_q), the last interval running from the last event time to the
#' horizon. The same sum as competing_risk's rsf_rmtl single-cause branch.
rsf_rmst <- function(o, horizon, oob) {
  if (o$family != "surv") {
    stop("rsf_rmst: expected a 'surv' family forest, got '", o$family, "'")
  }
  times <- o$time.interest
  n_t <- length(times)
  surv <- if (oob) o$survival.oob else o$survival
  surv <- matrix(surv, ncol = n_t)
  as.vector(horizon - (1 - surv) %*% c(diff(times), horizon - times[n_t]))
}

#' RMST at new rows from a fitted rsf_se_fit() forest
rsf_predict_rmst <- function(fit, X_new, horizon) {
  rsf_rmst(predict(fit, newdata = as.data.frame(X_new)), horizon, oob = FALSE)
}

#' The DR nuisance list, with the same fields as nuisance_pseudo_rf_oob's
rsf_se_nuisances <- function(pseudo, W, mu0, mu1, mu_cf, W.hat) {
  list(
    po = dr_pseudo(pseudo, W, mu1, mu0, W.hat),
    pseudo.hat = W * mu1 + (1 - W) * mu0,
    pseudo0.hat = mu0,
    pseudo.hat.cf = mu_cf,
    W.hat = W.hat
  )
}

#' Whole-sample OOB RSF nuisances ("oob")
#'
#' One forest per arm: each unit's own-arm prediction OOB, its other-arm
#' prediction from a forest that never saw it. pseudo.hat.cf from a pooled
#' forest on X alone (OOB).
nuisance_rsf_se_oob <- function(X, Y, D, W, horizon, pseudo, num.threads = NULL) {
  n_obs <- nrow(X)

  arm_fit <- function(arm) {
    in_arm <- W == arm
    fit <- rsf_se_fit(X[in_arm, , drop = FALSE], Y[in_arm], D[in_arm], horizon)
    pred <- numeric(n_obs)
    pred[in_arm] <- rsf_rmst(fit, horizon, oob = TRUE)
    pred[!in_arm] <- rsf_predict_rmst(fit, X[!in_arm, , drop = FALSE], horizon)
    pred
  }
  # control arm first, as t_learner_rf
  mu0 <- arm_fit(0)
  mu1 <- arm_fit(1)

  mu_cf <- rsf_rmst(rsf_se_fit(X, Y, D, horizon), horizon, oob = TRUE)

  W.hat <- trim_ps(predict(regression_forest(X, W, num.threads = num.threads))$predictions)

  rsf_se_nuisances(pseudo, W, mu0, mu1, mu_cf, W.hat)
}

#' Single leave-one-fold-out RSF nuisances ("scf")
#'
#' The crossfit twin of nuisance_rsf_se_oob: every forest is fit on the
#' training folds and predicts the held-out fold. Whole-sample pseudo-values
#' stay in the correction term, as in every DR arm.
nuisance_rsf_se_scf <- function(
  X,
  Y,
  D,
  W,
  horizon,
  pseudo,
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
        fit <- rsf_se_fit(X_train[rows, , drop = FALSE], Y[in_train][rows],
                          D[in_train][rows], horizon)
        rsf_predict_rmst(fit, X_test, horizon)
      }
      # control arm first, as nuisance_rsf_se_oob
      mu0 <- fit_predict(W_train == 0)
      mu1 <- fit_predict(W_train == 1)
      mu_cf <- fit_predict(rep(TRUE, sum(in_train)))

      W.hat.model <- regression_forest(X_train, W_train, num.threads = num.threads)
      W.hat <- trim_ps(predict(W.hat.model, newdata = X_test,
                               num.threads = num.threads)$predictions)

      c(list(fold = fold),
        rsf_se_nuisances(pseudo[in_test], W[in_test], mu0, mu1, mu_cf, W.hat))
    },
    .options = furrr_options(seed = TRUE)
  )

  scatter_folds(
    cross_fits,
    fold_indices,
    c("po", "pseudo.hat", "pseudo0.hat", "pseudo.hat.cf", "W.hat")
  )
}
