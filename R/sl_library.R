##########
# title: SuperLearner libraries, learner wrappers and the library pretest
##########
# Every SuperLearner call in the simulations draws its library from here, so the
# learners are chosen in one place. Depends only on SuperLearner and the learner
# packages (glmnet, gam, earth, ranger); sources nothing else from R/.
#
# ONE LIBRARY PER NUISANCE. The DR-learner fits three different regressions, and
# they used to share a single library (glm, glmnet, earth, gam, mean, ranger,
# with earth and ranger dropped at n <= 100):
#   W    the propensity. Every study is an RCT with P(W = 1) = 0.5, so flexible
#        learners only fit noise, which the 1 / (e(1 - e)) weight then
#        amplifies. mean + glm keeps the efficiency gain from estimating a known
#        propensity without the variance.
#   Y    the outcome model, fit SEPARATELY IN EACH ARM (a T-learner) on X alone,
#        so each fit sees about half the training rows: ~37 per arm at n = 100
#        (4 folds), ~100 at n = 250, ~225 at n = 500, ~450 at n = 1000. The
#        library is sized to that: at n <= 100 no glm or gam (11+ parameters on
#        ~37 rows, ~15 events for a binary outcome), glm and gam from n = 250,
#        earth from n = 500.
#   tau  the stage-2 regression of the pseudo-outcome on X. Chosen without
#        reference to the DGMs' HTE shapes: learners from each function class
#        (constant, linear, penalised linear, penalised pairwise interactions,
#        additive smooth, adaptive splines, tree ensemble), regularised because
#        the pseudo-outcome is noisy - in an RCT its residual is doubled, so its
#        variance is about 4x. The interaction lasso and earth join at n >= 250.
#        Each lasso comes at BOTH lambda.min and lambda.1se, and SuperLearner's
#        CV weighs them: fixing either in advance is a tuning choice, and
#        making it from how the scenarios came out would be tuning to the DGM.
#        (lambda.1se alone shrank the strong linear effect of scenario 2 and
#        cost ~3-4x in MSE at n >= 500; with both, CV put 55-67% of the weight
#        on lambda.min there. Where there is no heterogeneity lambda.1se
#        shrinks to a constant, i.e. to SL.mean, so it costs nothing.) A
#        10-learner version with ridge, a second gam and BART did worse than
#        the old library at n = 500, from noisier ensemble weights.
#
# CUSTOM WRAPPERS AND WORKERS. SuperLearner looks learners up by name with
# get(name, envir = env). The studies fit inside future_map() on multisession
# workers, where a wrapper defined in the parent's global environment does not
# exist and a name given as a string is not detected as a global. So the
# wrappers live in SL_ENV, which every call here passes as `env =`; code that
# references SL_ENV by symbol makes future export it. Its parent is the
# SuperLearner namespace, so the built-in SL.* and screen.All resolve through
# it too.
#
# PREDICTION. SL.glmnet.int and SL.glmnet.int.1se expand X into interactions,
# so predict() on them would be handed the unexpanded columns. Fit with
# sl_fit_predict(), which predicts through SuperLearner's own newX at fit time;
# their predict methods stop() rather than mispredict.

require(SuperLearner)

# ---- libraries --------------------------------------------------------------

#' Per-nuisance SuperLearner libraries for a study of size n
#'
#' @param n total sample size of the dataset the libraries will be fit to
#' @return list(W = , Y = , tau = ) of learner names
sl_libraries <- function(n) {
  list(
    W = c("SL.mean", "SL.glm"),
    Y = c("SL.mean", "SL.glmnet", "SL.ranger.ns25",
          if (n > 100) c("SL.glm", "SL.gam"),
          if (n >= 500) "SL.earth"),
    tau = c("SL.mean", "SL.glm", "SL.glmnet", "SL.glmnet.1se", "SL.gam",
            "SL.ranger.ns25",
            if (n > 100) c("SL.glmnet.int", "SL.glmnet.int.1se", "SL.earth"))
  )
}

#' Normalise an sl_lib argument to list(W = , Y = , tau = )
#'
#' A plain character vector is the old single-library form and is used for
#' all three nuisances, so callers that still pass one (R/regression_check.R's
#' SL_LIB, crossfitting/confidence_intervals/cf_ci_testing.R) behave as before.
as_sl_libs <- function(sl_lib) {
  if (is.null(sl_lib)) return(NULL)
  if (is.character(sl_lib)) return(list(W = sl_lib, Y = sl_lib, tau = sl_lib))
  missing_parts <- setdiff(c("W", "Y", "tau"), names(sl_lib))
  if (length(missing_parts)) {
    stop("sl_lib is missing ", paste(missing_parts, collapse = ", "),
         " - pass sl_libraries(n) or a character vector")
  }
  sl_lib
}

# ---- learner wrappers -------------------------------------------------------

SL_ENV <- local({

  # the lasso at lambda.1se; SL.glmnet itself is lambda.min (useMin = TRUE).
  # Both go in the stage-2 library - see the header.
  SL.glmnet.1se <- function(Y, X, newX, family, obsWeights, ...) {
    SuperLearner::SL.glmnet(Y, X, newX, family, obsWeights, useMin = FALSE, ...)
  }

  # lasso over every main effect and pairwise interaction, at lambda.min
  # (SL.glmnet.int) or lambda.1se (SL.glmnet.int.1se)
  glmnet_int <- function(Y, X, newX, family, obsWeights, useMin, cls, ...) {
    f <- ~ .^2
    out <- SuperLearner::SL.glmnet(
      Y, stats::model.matrix(f, X)[, -1, drop = FALSE],
      stats::model.matrix(f, newX)[, -1, drop = FALSE],
      family, obsWeights, useMin = useMin, ...
    )
    class(out$fit) <- cls
    out
  }
  SL.glmnet.int <- function(Y, X, newX, family, obsWeights, ...) {
    glmnet_int(Y, X, newX, family, obsWeights, useMin = TRUE, cls = "SL.glmnet.int", ...)
  }
  SL.glmnet.int.1se <- function(Y, X, newX, family, obsWeights, ...) {
    glmnet_int(Y, X, newX, family, obsWeights, useMin = FALSE,
               cls = "SL.glmnet.int.1se", ...)
  }

  # ranger's regression default, min.node.size = 5, overfits the pseudo-outcome
  SL.ranger.ns25 <- function(Y, X, newX, family, obsWeights, ...) {
    SuperLearner::SL.ranger(Y, X, newX, family, obsWeights, min.node.size = 25, ...)
  }

  environment()
}, envir = new.env(parent = asNamespace("SuperLearner")))

no_predict <- function(object, ...) {
  stop(class(object)[1], " cannot predict after fitting - fit it with ",
       "sl_fit_predict(), which predicts at newX (see R/sl_library.R)")
}
predict.SL.glmnet.int <- no_predict
predict.SL.glmnet.int.1se <- no_predict

# ---- fitting ----------------------------------------------------------------

#' Fit a SuperLearner and predict at one or more covariate sets
#'
#' Predicts through SuperLearner's own newX, at fit time, so no learner needs
#' a predict method: the newX_list pieces are stacked, predicted in one fit,
#' and split back. Gives the same predictions as fitting and then calling
#' predict(fit, newdata = ) for every built-in learner.
#'
#' FAILED FITS. If SuperLearner() itself errors, every row is predicted as the
#' (weighted) mean of Y, with a warning, and `failed` carries the error. The
#' pretest's SL.mean fallback (bug K) cannot cover this: a learner can pass the
#' pretest's 2-fold CV and still fail the live 10-fold one. It happens with
#' per-arm binary outcome models at n = 100, where an arm's training rows can
#' hold ~2 events - binary scenario 4, whose treated risk is low - so every
#' learner errors on some inner fold ("All algorithms dropped from library").
#'
#' ZERO WEIGHTS. When NNLS gives every learner weight 0, SuperLearner only
#' warns ("All algorithms have zero weight") and SL.predict is 0 on every row -
#' not an estimate of anything. It happens to the DR-learner's stage 2 on a
#' binary outcome: the risk-difference effect is small next to the
#' pseudo-outcome's noise, so every learner's CV predictions, SL.mean's
#' leave-fold-out mean above all, are uncorrelated or negatively correlated
#' with the held-out pseudo-outcomes, and NNLS has no intercept to fall back
#' on. The fit then uses the discrete SuperLearner - the single learner with
#' the lowest CV risk - with a warning, and `zero_weights` names it.
#'
#' @param newX_list named list of covariate data.frames, each with X's columns
#' @return list(pred = named list of numeric vectors, coef = ensemble weights
#'   (NULL after a failed fit), failed = NA or the error message,
#'   zero_weights = NA or the learner used because every weight was 0)
sl_fit_predict <- function(Y, X, newX_list, SL.library, family = gaussian(),
                           obsWeights = NULL, cvControl = list()) {
  sizes <- vapply(newX_list, NROW, integer(1))
  newX <- do.call(rbind, unname(newX_list))
  method <- if (identical(family$family, "binomial")) "method.NNloglik" else "method.NNLS"
  piece <- factor(rep(names(newX_list), sizes), levels = names(newX_list))

  fit <- tryCatch(
    SuperLearner(Y = Y, X = X, newX = newX, family = family,
                 SL.library = as.character(SL.library), method = method,
                 obsWeights = obsWeights, cvControl = cvControl, env = SL_ENV),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    msg <- conditionMessage(fit)
    warning("sl_fit_predict: SuperLearner failed (", msg, "); predicting the mean of Y.")
    w <- if (is.null(obsWeights)) rep(1, length(Y)) else obsWeights
    pred <- rep(stats::weighted.mean(Y, w), sum(sizes))
    return(list(pred = split(pred, piece), coef = NULL, failed = msg,
                zero_weights = NA_character_))
  }

  pred <- as.numeric(fit$SL.predict)
  zero_weights <- NA_character_
  if (all(fit$coef == 0)) {
    ok <- is.finite(fit$cvRisk)
    if (length(fit$errorsInLibrary)) ok <- ok & !fit$errorsInLibrary
    if (any(ok)) {
      best <- which(ok)[which.min(fit$cvRisk[ok])]
      zero_weights <- names(fit$cvRisk)[best]
      pred <- as.numeric(fit$library.predict[, best])
    } else {
      zero_weights <- "(none usable, used the mean)"
      w <- if (is.null(obsWeights)) rep(1, length(Y)) else obsWeights
      pred <- rep(stats::weighted.mean(Y, w), sum(sizes))
    }
    warning("sl_fit_predict: every ensemble weight is zero; using ", zero_weights, ".")
  }
  list(pred = split(pred, piece), coef = fit$coef, failed = NA_character_,
       zero_weights = zero_weights)
}

#' Record a failed or zero-weight sl_fit_predict() fit in a pretested library's
#' "dropped" attribute, so dropped_table() reports it alongside the pretest's
#' drops
mark_failed_fit <- function(lib, fit) {
  if (!is.na(fit$failed)) {
    attr(lib, "dropped") <- c(attr(lib, "dropped"),
                              "(whole fit)" = paste("SuperLearner failed, used the mean:",
                                                    fit$failed))
  }
  if (!is.null(fit$zero_weights) && !is.na(fit$zero_weights)) {
    attr(lib, "dropped") <- c(attr(lib, "dropped"),
                              "(whole fit)" = paste("all ensemble weights zero, used",
                                                    fit$zero_weights))
  }
  lib
}

# ---- library pretest --------------------------------------------------------

#' Drop SuperLearner algorithms that error or give non-finite predictions
#'
#' Fits each candidate on its own with a 2-fold inner CV and keeps the
#' survivors. Warnings are recorded, not fatal: this used to catch them with
#' tryCatch(warning = ), which aborted the fit at the first warning and dropped
#' the learner, so a benign glm or gam warning emptied a fold's library down to
#' one or two learners (see bug L in R/cate_models.R).
#'
#' @return the surviving library, with attr "dropped" (named reasons for the
#'   removed learners) and attr "warned" (the first warning of each kept
#'   learner that warned)
pretest_superlearner <- function(Y, X, SL.library, family) {
  working_lib <- character()
  dropped <- character()
  warned <- character()
  for (alg in SL.library) {
    warnings_seen <- character()
    fit <- tryCatch(
      withCallingHandlers(
        SuperLearner(Y = Y, X = X, SL.library = alg, family = family,
                     cvControl = list(V = 2), env = SL_ENV),
        warning = function(w) {
          warnings_seen <<- c(warnings_seen, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )
    reason <- if (inherits(fit, "error")) {
      paste("error:", conditionMessage(fit))
    } else if (any(fit$errorsInCVLibrary) || any(fit$errorsInLibrary)) {
      "error inside SuperLearner"
    } else if (is.null(fit$SL.predict) || !all(is.finite(fit$SL.predict))) {
      "non-finite predictions"
    } else {
      NA_character_
    }
    if (is.na(reason)) {
      working_lib <- c(working_lib, alg)
      if (length(warnings_seen)) warned[alg] <- warnings_seen[1]
    } else {
      dropped[alg] <- reason
    }
  }
  if (length(dropped) > 0) {
    cat("Removed libraries due to NA/error:\n")
    print(dropped)
  }
  if (length(working_lib) == 0) {
    # bug K: every candidate failed on this fold - falling through with
    # character(0) sends an empty SL.library into the caller's live SuperLearner()
    # call, which crashes building a 0-column predictions data.frame() ("arguments
    # imply differing number of rows: 0, 1"). SL.mean is asserted directly, not
    # re-run through the loop above, so it can't recursively trigger the same
    # emptying it's meant to prevent.
    warning("pretest_superlearner: every candidate failed on this fold; ",
            "falling back to SL.mean.")
    working_lib <- "SL.mean"
  }
  structure(working_lib, dropped = dropped, warned = warned)
}

#' Dropped learners from a set of pretested libraries, as one data.frame
#'
#' @param libs named list of pretest_superlearner() results, e.g. list(Y = , W = )
#' @param fold the fold they were pretested on
dropped_table <- function(libs, fold) {
  rows <- lapply(names(libs), function(model) {
    d <- attr(libs[[model]], "dropped")
    if (!length(d)) return(NULL)
    data.frame(fold = fold, model = model, learner = names(d), reason = unname(d),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
