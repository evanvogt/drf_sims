##########
# title: CATE model fitting functions - missing covariates and continuous outcome
##########

# packages

require(grf)
require(SuperLearner)
require(dplyr)
require(future)
library(furrr)

# Main function
run_all_cate_methods <- function(data, n_folds = 10, workers = 2, sl_lib = NULL, fmla_info = NULL, ipw = NULL) {
  
  # Set up parallelisation
  metaplan <- plan(multisession, workers = workers)
  on.exit(plan(metaplan), add = TRUE)
  
  # Data set up
  X <- as.matrix(data[, -c(1:2)])
  Y <- data$Y
  W <- data$W
  
  # fold creation
  n_obs <- nrow(X)
  fold_indices <- sort(seq(n_obs) %% n_folds) + 1
  fold_list <- unique(fold_indices)
  fold_pairs <- utils::combn(fold_list, 2, simplify = F)
  
  # results container
  results <- list()
  
  # nuisance estimation (RF)
  cat("Computing nuisance functions...\n")
  nuisances_rf <- nuisance_rf(X, Y, W, fold_indices, fold_pairs, ipw)
  
  # causal forest
  cat("Running Causal Forest...\n")
  results$causal_forest <- run_causal_forest(X, Y, W, nuisances_rf, fold_indices, fold_list, ipw)
  
  # DR with Random forest (separate function not required bc we have stage 1 and 2 functions)
  cat("Running DR Random Forest...\n")
  results$dr_random_forest <- list(tau = stage_2_rf(X, nuisances_rf$po_matrix, fold_indices, fold_list, ipw))
  
  
  if(!anyNA(X)) {
    # DR oracle
    if (!is.null(fmla_info)) {
      cat("Running DR Oracle...\n")
      
      results$dr_oracle <- run_dr_oracle(X, Y, W, fmla_info, fold_indices, fold_list, ipw)
    }
    
    # DR semi oracle
    cat("Running DR Semi-Oracle...\n")
    results$dr_semi_oracle <- run_dr_semi_oracle(X, Y, W, fold_indices, fold_list, ipw)
    
    # DR SuperLearner (if data is complete)
    cat("Running DR Super Learner...")
    X <- as.data.frame(X)
    nuisances_sl <- nuisance_sl(X, Y, W, fold_indices, fold_pairs, sl_lib, ipw)
    
    results$dr_superlearner <- list(tau = stage_2_sl(X, nuisances_sl$po_matrix, fold_indices, fold_list, sl_lib, ipw))
    results$nuisances_sl <- nuisances_sl
  }
  
  # keep nuisances for running post estimation tests
  results$nuisances_rf <- nuisances_rf
  results$fold_indices <- fold_indices
  
  return(results)
}

# fit the causal forest
run_causal_forest <- function(X, Y, W, nuisances, fold_indices, fold_list, ipw = NULL) {
  n_obs <- nrow(X)
  
  # Get tau predictions using causal forest with pre-computed nuisance functions
  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold <- !in_train
    
    # Use nuisance predictions that didn't use the current fold
    forest <- causal_forest(X[in_train,], Y[in_train], W[in_train], 
                            nuisances$Y.hat.cf_matrix[in_train, fold], 
                            nuisances$W.hat_matrix[in_train, fold],
                            sample.weights = if(!is.null(ipw)) ipw[in_train] else NULL)
    
    pred <- predict(forest, newdata = X[in_fold,], estimate.variance = TRUE)
    
    list(fold = fold, tau = pred$predictions, variance = pred$variance.estimates)
  }, .options = furrr_options(seed = TRUE))
  
  # Reconstruct tau vector
  tau <- rep(NA, n_obs)
  tau_var <- rep(NA, n_obs)
  
  for (result in tau_results) {
    fold <- result$fold
    in_fold <- fold_indices == fold
    tau[in_fold] <- result$tau
    tau_var[in_fold] <- result$variance
  }
  
  return(list(
    tau = tau,
    variance = tau_var
  ))
}

# fit the DR Oracle
run_dr_oracle <- function(X, Y, W, fmla_info, fold_indices, fold_list, ipw = NULL) {
  n_obs <- nrow(X)
  
  # Get parameters for oracle
  X <- as.data.frame(X)
  list2env(fmla_info$params, envir = environment())
  fmla <- parse(text = fmla_info$fmla)
  
  # calculate true values with formula
  W_temp <- rep(1, n_obs)
  Y1.hat <- eval(fmla, envir = list2env(c(list(W = W_temp), X)))
  
  W_temp <- rep(0, n_obs)
  Y0.hat <- eval(fmla, envir = list2env(c(list(W = W_temp), X)))
  
  Y.hat <- eval(fmla, envir = list2env(c(list(W = W), X)))
  W.hat <- rep(0.5, n_obs)
  
  X <- as.matrix(X)
  
  # Calculate pseudo-outcomes
  cate <- Y1.hat - Y0.hat
  po <- cate + ((Y - Y.hat) * (W - W.hat)) / (W.hat * (1 - W.hat))
  
  tau <- stage_2_rf(X, po, fold_indices, fold_list, ipw)
  return(list(
    tau = tau
  ))
}

# fit DR semi-oracle (fixed W.hat)
run_dr_semi_oracle <- function(X, Y, W, fold_indices, fold_list, ipw = NULL) {
  n_obs <- nrow(X)
  
  # set W.hat
  W.hat <- rep(0.5, n_obs)
  
  
  # Parallel cross-fitting for Y model only
  cross_fits <- future_map(seq_along(fold_list), function(fold) {
    in_train <- fold_indices != fold
    in_test <- !in_train
    
    Y.hat.model <- regression_forest(cbind(W[in_train], X[in_train, ]), Y[in_train], sample.weights = if(!is.null(ipw)) ipw[in_train] else NULL)
    
    X_test <- X[in_test, ]
    Y0.hat <- predict(Y.hat.model, newdata = cbind(W = 0, X_test))$predictions
    Y1.hat <- predict(Y.hat.model, newdata = cbind(W = 1, X_test))$predictions
    
    W_test <- W[in_test]
    Y.hat <- W_test * Y1.hat + (1 - W_test) * Y0.hat
    cate <- Y1.hat - Y0.hat
    po <- cate + ((Y[in_test] - Y.hat) * (W_test - 0.5)) / (0.5 * 0.5)
    
    list(po = po, Y0.hat = Y0.hat, fold = fold)
  }, .options = furrr_options(seed = TRUE))
  
  # Reorganize results
  po <- rep(NA, n_obs)
  Y0.hat <- rep(NA, n_obs)
  
  for (i in seq_along(cross_fits)) {
    fold <- cross_fits[[i]]$fold
    in_fold <- fold_indices == fold
    po[in_fold] <- cross_fits[[i]]$po
    Y0.hat[in_fold] <- cross_fits[[i]]$Y0.hat
  }
  
  tau <- stage_2_rf(X, po, fold_indices, fold_list, ipw)
  return(list(
    tau = tau,
    po = po,
    Y0.hat = Y0.hat
  ))
}

# second stage regression function (RF)
stage_2_rf <- function(X, po, fold_indices, fold_list, ipw = NULL) {
  n_obs <- nrow(X)
  
  # check if we have a po_matrix or vector
  single <- is.vector(po)
  
  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold <- !in_train
    
    if (single) {
      # Single vector case (oracle, semi-oracle)
      forest <- regression_forest(X[in_train, ], po[in_train], sample.weights = if(!is.null(ipw)) ipw[in_train] else NULL)
    } else {
      # Matrix case (RF)
      forest <- regression_forest(X[in_train, ], po[in_train, fold], sample.weights = if(!is.null(ipw)) ipw[in_train] else NULL)
    }
    tau_pred <- predict(forest, newdata = X[in_fold, ])$predictions
    list(fold = fold, predictions = tau_pred)
  }, .options = furrr_options(seed = TRUE))
  
  # Reconstruct tau vector
  tau <- rep(NA, n_obs)
  for (result in tau_results) {
    fold <- result$fold
    in_fold <- fold_indices == fold
    tau[in_fold] <- result$predictions
  }
  return(tau)
}

# second stage regression function (SuperLearner)
stage_2_sl <- function(X, po, fold_indices, fold_list, sl_lib, ipw = NULL) {
  n_obs <- nrow(X)
  
  # check if we have a po_matrix or vector
  single <- is.vector(po)
  
  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold <- !in_train
    
    if (single) {
      po_lib <- pretest_superlearner(po[in_train], X[in_train,], sl_lib, gaussian())
      po_model <- SuperLearner(po[in_train], X[in_train, ], family = gaussian(), SL.library = po_lib, obsWeights = if(!is.null(ipw)) ipw[in_train] else NULL)
    } else {
      po_lib <- pretest_superlearner(po[in_train, fold], X[in_train,], sl_lib, gaussian())
      po_model <- SuperLearner(po[in_train, fold], X[in_train, ], family = gaussian(), SL.library = sl_lib, obsWeights = if(!is.null(ipw)) ipw[in_train] else NULL)
    }
    
    tau_pred <- predict(po_model, newdata = X[in_fold, ])$pred
    list(fold = fold, predictions = tau_pred)
  }, .options = furrr_options(seed = TRUE))
  
  # Reconstruct tau vector
  tau <- rep(NA, n_obs)
  for (result in tau_results) {
    fold <- result$fold
    in_fold <- fold_indices == fold
    tau[in_fold] <- result$predictions
  }
  return(tau)
}

# crossing fitting of nuisance parameters using RF
nuisance_rf <- function(X, Y, W, fold_indices, fold_pairs, ipw = NULL) {
  
  # double crossfitting across fold pairs
  cross_fits <- future_map(seq_along(fold_pairs), function(i) {
    fold_pair <- fold_pairs[[i]]
    in_train <- !(fold_indices %in% fold_pair)
    in_test <- !in_train
    
    Y.hat.model <- regression_forest(cbind(W[in_train], X[in_train, ]), Y[in_train], sample.weights = if(!is.null(ipw)) ipw[in_train] else NULL)
    Y.hat.cf.model <- regression_forest(X[in_train, ], Y[in_train], sample.weights = if(!is.null(ipw)) ipw[in_train] else NULL)
    W.hat.model <- regression_forest(X[in_train, ], W[in_train], sample.weights = if(!is.null(ipw)) ipw[in_train] else NULL)
    
    X_test <- X[in_test, ]
    
    Y0.hat <- predict(Y.hat.model, newdata = cbind(W = 0, X_test))$predictions
    Y1.hat <- predict(Y.hat.model, newdata = cbind(W = 1, X_test))$predictions
    Y.hat.cf <- predict(Y.hat.cf.model, newdata = X_test)$predictions
    W.hat <- predict(W.hat.model, newdata = X_test)$predictions
    
    W_test <- W[in_test]
    Y.hat <- W_test * Y1.hat + (1 - W_test) * Y0.hat
    
    cate <- Y1.hat - Y0.hat
    po <- cate + ((Y[in_test] - Y.hat) * (W_test - W.hat)) / (W.hat * (1 - W.hat))
    
    list(po = po, Y.hat = Y.hat, Y0.hat = Y0.hat, Y.hat.cf = Y.hat.cf, W.hat = W.hat, fold_pair = fold_pair)
  }, .options = furrr_options(seed = TRUE))
  
  # make matrices of predctions
  fold_list <- unique(fold_indices)
  po_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "po")
  Y.hat_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "Y.hat")
  Y.hat.cf_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "Y.hat.cf")
  Y0.hat_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "Y0.hat") 
  W.hat_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "W.hat")
  
  return(list(
    po_matrix = po_matrix,
    Y.hat_matrix = Y.hat_matrix,
    Y.hat.cf_matrix = Y.hat.cf_matrix,
    Y0.hat_matrix = Y0.hat_matrix,
    W.hat_matrix = W.hat_matrix
  ))
}

# crossfitting for nuisances with SL
nuisance_sl <- function(X, Y, W, fold_indices, fold_pairs, sl_lib, ipw = NULL) {
  
  # double crossfitting across fold pairs
  cross_fits <- future_map(seq_along(fold_pairs), function(i) {
    fold_pair <- fold_pairs[[i]]
    in_train <- !(fold_indices %in% fold_pair)
    in_test <- !in_train
    
    X_train <- X[in_train, ]
    X_W_train <- cbind(W = W[in_train], X_train)
    
    # train models + pretest SL libraries
    Y_lib <- pretest_superlearner(Y[in_train], X_W_train, sl_lib, gaussian())
    Y.hat.model <- SuperLearner(Y = Y[in_train], X = X_W_train, SL.library = Y_lib, obsWeights = if(!is.null(ipw)) ipw[in_train] else NULL)
    
    W_lib <- pretest_superlearner(W[in_train], X_train, sl_lib, binomial())
    W.hat.model <- SuperLearner(W[in_train], X_train, family = binomial(), 
                                SL.library = W_lib, method = method.NNloglik(), obsWeights = if(!is.null(ipw)) ipw[in_train] else NULL)
    
    # predictions
    X_test <- X[in_test, ]
    Y0.hat <- predict(Y.hat.model, newdata = cbind(W = 0, X_test))$pred
    Y1.hat <- predict(Y.hat.model, newdata = cbind(W = 1, X_test))$pred
    W.hat <- predict(W.hat.model, newdata = X_test)$pred
    
    # failsafe is the SuperLearner fails
    if (all(Y0.hat == 0) && all(Y1.hat == 0)) {
      warning("SuperLearner failed for Y.hat. Using mean(Y).")
      Y0.hat <- rep(mean(Y[in_train][W[in_train] == 0], na.rm = TRUE), sum(in_test))
      Y1.hat <- rep(mean(Y[in_train][W[in_train] == 1], na.rm = TRUE), sum(in_test))
    }
    
    if (all(W.hat == 0)) {
      warning("SuperLearner failed for W.hat. Using mean(W).")
      W.hat <- rep(mean(W[in_train], na.rm = TRUE), sum(in_test))
    }
    
    # Trim extreme propensity scores
    W.hat[W.hat < 0.05] <- 0.05
    W.hat[W.hat > 0.95] <- 0.95
    
    # pseudo outcomes
    W_test <- W[in_test]
    Y.hat <- W_test * Y1.hat + (1 - W_test) * Y0.hat
    
    cate <- Y1.hat - Y0.hat
    po <- cate + ((Y[in_test] - Y.hat) * (W_test - W.hat)) / (W.hat * (1 - W.hat))
    
    list(po = po, Y.hat = Y.hat, Y0.hat = Y0.hat, W.hat = W.hat, fold_pair = fold_pair)
  }, .options = furrr_options(seed = TRUE))
  
    # make matrices of predctions
  fold_list <- unique(fold_indices)
  po_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "po")
  Y.hat_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "Y.hat")
  Y0.hat_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "Y0.hat") 
  W.hat_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "W.hat")
  
  return(list(
    po_matrix = po_matrix,
    Y.hat_matrix = Y.hat_matrix,
    Y0.hat_matrix = Y0.hat_matrix,
    W.hat_matrix = W.hat_matrix
  ))
}

# collation of list to matrix
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

# check superlearner algorithms and remove failures
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
  return(working_lib)
}