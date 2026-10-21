##########
# title: SuperLearner libraries, learner wrappers and the library pretest
##########
# Every SuperLearner call in the simulations takes its library from here.
# Depends only on SuperLearner and the learner packages; sources nothing from R/.
#
# ONE LIBRARY PER NUISANCE:
#   W    propensity: mean + glm. Every study is an RCT at P(W = 1) = 0.5, so a
#        flexible learner only fits noise that the 1 / (e(1 - e)) weight amplifies.
#   Y    outcome, fit per arm (T-learner), so each fit sees about half the
#        training rows (~37 at n = 100). glm and gam join at n = 250, earth at 500.
#   tau  stage 2: regularised learners from each function class, chosen without
#        reference to the DGMs (the pseudo-outcome's variance is ~4x). Each
#        lasso comes at both lambda.min and lambda.1se and CV weighs them, since
#        fixing one would be tuning to the DGM (lambda.1se alone cost ~3-4x MSE
#        in scenario 2 at n >= 500).
#
# WORKERS. SuperLearner finds learners with get(name, envir = env), and a
# wrapper in the parent's global environment does not exist on a multisession
# worker. So the wrappers live in SL_ENV (parent: the SuperLearner namespace),
# passed as `env =` everywhere here.
#
# PREDICTION. SL.glmnet.int(.1se) expand X into interactions, so fit them with
# sl_fit_predict(), which predicts via newX at fit time; their predict methods stop().

library(SuperLearner)

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
#' A plain character vector is used for all three nuisances (R/regression_check.R's
#' SL_LIB and crossfitting/confidence_intervals/cf_ci_testing.R pass one).
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
#' pretest's SL.mean fallback cannot cover this: a learner can pass the
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
#' survivors. Warnings are recorded, not fatal: tryCatch(warning = ) would
#' abort the fit at the first warning and drop the learner, so a benign glm or
#' gam warning would empty a fold's library down to one or two learners.
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
    # every candidate failed on this fold - falling through with
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
