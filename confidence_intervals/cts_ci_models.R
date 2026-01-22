##########
# title: CATE model fitting with confidence intervals - cts outcome
##########

# packages

require(grf)
require(SuperLearner)
require(dplyr)
require(future)
require(furrr)

# Main function
run_all_cate_methods <- function(data, n_folds = 10, fmla_info = NULL, CI_boot = 200, CI_sf = 0.5, alpha = 0.05) {
  
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
  nuisances_rf <- nuisance_rf(X, Y, W, fold_indices, fold_pairs)
  
  # causal forest
  cat("Running Causal Forest...\n")
  results$causal_forest <- run_causal_forest(X, Y, W, nuisances_rf, fold_indices, fold_list)
  
  # bootstrap CIs
  cat("Running Causal Forest bootstrap... \n")
  results$causal_forest <- c(results$causal_forest, cf_half_boot(X, Y, W, nuisances_rf, results$causal_forest$tau, CI_boot, CI_sf, alpha, fold_indices, fold_list))
  
  # DR with Random forest (separate function not required bc we have stage 1 and 2 functions)
  cat("Running DR Random Forest...\n")
  results$dr_random_forest <- list(tau = stage_2_rf(X, nuisances_rf$po_matrix, fold_indices, fold_list))
  
  cat("Running DR RF bootstrap... \n")
  results$dr_random_forest <- c(results$dr_random_forest, rf_half_boot(X, Y, W, nuisances_rf$po_matrix, results$dr_random_forest$tau, CI_boot, CI_sf, alpha, fold_indices, fold_list))
  
  
  # DR oracle
  if (!is.null(fmla_info)) {
    cat("Running DR Oracle...\n")
    
    results$dr_oracle <- run_dr_oracle(X, Y, W, fmla_info, fold_indices, fold_list)
    
    cat("Runnings Oracle bootstrap...\n")
    results$dr_oracle <- c(results$dr_oracle, rf_half_boot(X, Y, W, results$dr_oracle$po, results$dr_oracle$tau, CI_boot, CI_sf, alpha, fold_indices, fold_list))
  }
  
  # DR semi oracle
  cat("Running DR Semi-Oracle...\n")
  results$dr_semi_oracle <- run_dr_semi_oracle(X, Y, W, fold_indices, fold_list)
  
  cat("Running Semi-Oracle bootstrap...\n")
  results$dr_semi_oracle <- c(results$dr_semi_oracle, rf_half_boot(X, Y, W, results$dr_semi_oracle$po, results$dr_semi_oracle$tau, CI_boot, CI_sf, alpha, fold_indices, fold_list))
  
  # keep nuisances for running post estimation tests
  results$nuisances_rf <- nuisances_rf
  results$fold_indices <- fold_indices
  
  return(results)
}

# fit the causal forest - use in built CI generation
run_causal_forest <- function(X, Y, W, nuisances, fold_indices, fold_list) {
  n_obs <- nrow(X)
  
  # Get tau predictions using causal forest with pre-computed nuisance functions
  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold <- !in_train
    
    # Use nuisance predictions that didn't use the current fold
    forest <- causal_forest(X[in_train,], Y[in_train], W[in_train], 
                            nuisances$Y.hat.cf_matrix[in_train, fold], 
                            nuisances$W.hat_matrix[in_train, fold])
    
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

# half-sample bootstrap for CF
cf_half_boot <- function(X, Y, W, nuisances, tau, CI_boot = 200, CI_sf = 0.5, alpha = 0.05, fold_indices, fold_list) {
  
  n_obs <- nrow(X)
  fold_membership <- lapply(fold_list, function(fold) which(fold_indices == fold))
  fold_sizes <- lengths(fold_membership)
  
  draws <- future_map(seq_len(CI_boot), function(b) {
    # half sample construction
    half_samples <- rep(F, n_obs)
    for (fold in fold_list) {
      fold_obs <- fold_membership[[fold]]
      n_half <- floor(fold_sizes[fold]/2)
      
      selected <- sample(fold_obs, n_half, replace = F)
      half_samples[selected] <- TRUE
    }
    
    # get half sample CATE
    tau_half_results <- lapply(fold_list, function(fold) {
      in_train <- half_samples & (fold_indices != fold)
      in_fold <- fold_indices == fold
      
      half_cf <- causal_forest(X[in_train,], Y[in_train], W[in_train], 
                               nuisances$Y.hat.cf_matrix[in_train, fold], 
                               nuisances$W.hat_matrix[in_train, fold],
                               sample.fraction = CI_sf)
      tau_half_est <- predict(half_cf, newdata = X[in_fold,])$predictions
      return(tau_half_est)
    })
    
    tau_half <- rep(NA, nrow(X))
    for (i in fold_list) {
      tau_half[fold_membership[[i]]] <- tau_half_results[[i]]
    }
    
    
    # construct root
    half_root <- tau - tau_half
    return(half_root)
  }, .options = furrr_options(seed = T))
  
  draws <- do.call(cbind, draws)
  
  # get quantiles of the half root distribution
  lambda_hat <- apply(draws, 1, var)
  draws_norm <- abs(draws)/sqrt(lambda_hat)
  col_max <- apply(draws_norm, 2, max)
  S_star <- quantile(col_max, 1 - (alpha/2))
  
  margin <- sqrt(lambda_hat)*S_star
  
  return(list(
    hb_lb = tau - margin,
    hb_ub = tau + margin,
    draws = draws
  ))
}

# half-sample bootstrap for DR-RF
rf_half_boot <- function(X, Y, W, po, tau, CI_boot = 200, CI_sf = 0.5, alpha = 0.05, fold_indices, fold_list) {
  
  n_obs <- nrow(X)
  fold_membership <- lapply(fold_list, function(fold) which(fold_indices == fold))
  fold_sizes <- lengths(fold_membership)
  
  single <- is.vector(po)
  
  draws <- future_map(seq_len(CI_boot), function(b) {
    # half sample construction
    half_samples <- rep(F, n_obs)
    for (fold in fold_list) {
      fold_obs <- fold_membership[[fold]]
      n_half <- floor(fold_sizes[fold]/2)
      
      selected <- sample(fold_obs, n_half, replace = F)
      half_samples[selected] <- TRUE
    }
    
    # get half sample CATE
    tau_half_results <- lapply(fold_list, function(fold) {
      in_train <- half_samples & (fold_indices != fold)
      in_fold <- fold_indices == fold
      
      if (single) {
        # Single vector case (oracle, semi-oracle)
        half_rf <- regression_forest(X[in_train, ], po[in_train], sample.fraction = CI_sf)
      } else {
        # Matrix case (RF)
        half_rf <- regression_forest(X[in_train, ], po[in_train, fold], sample.fraction = CI_sf)
      }

      tau_half_est <- predict(half_rf, newdata = X[in_fold,])$predictions
      return(tau_half_est)
    })
    
    tau_half <- rep(NA, nrow(X))
    for (i in fold_list) {
      tau_half[fold_membership[[i]]] <- tau_half_results[[i]]
    }
    
    
    # construct root
    half_root <- tau - tau_half
    return(half_root)
  }, .options = furrr_options(seed = T))
  
  draws <- do.call(cbind, draws)
  
  # get quantiles of the half root distribution
  lambda_hat <- apply(draws, 1, var)
  draws_norm <- abs(draws)/sqrt(lambda_hat)
  col_max <- apply(draws_norm, 2, max)
  S_star <- quantile(col_max, 1 - (alpha/2))
  
  margin <- sqrt(lambda_hat)*S_star
  
  return(list(
    hb_lb = tau - margin,
    hb_ub = tau + margin,
    draws = draws
  ))
}
# fit the DR Oracle
run_dr_oracle <- function(X, Y, W, fmla_info, fold_indices, fold_list) {
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
  
  tau <- stage_2_rf(X, po, fold_indices, fold_list)
  return(list(
    tau = tau,
    po = po,
    Y0.hat = Y0.hat
  ))
}

# fit DR semi-oracle (fixed W.hat)
run_dr_semi_oracle <- function(X, Y, W, fold_indices, fold_list) {
  n_obs <- nrow(X)
  
  # set W.hat
  W.hat <- rep(0.5, n_obs)
  
  
  # Parallel cross-fitting for Y model only
  cross_fits <- future_map(seq_along(fold_list), function(fold) {
    in_train <- fold_indices != fold
    in_test <- !in_train
    
    Y.hat.model <- regression_forest(cbind(W[in_train], X[in_train, ]), Y[in_train])
    
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
  
  tau <- stage_2_rf(X, po, fold_indices, fold_list)
  return(list(
    tau = tau,
    po = po,
    Y0.hat = Y0.hat
  ))
}

# second stage regression function (RF)
stage_2_rf <- function(X, po, fold_indices, fold_list) {
  n_obs <- nrow(X)
  
  # check if we have a po_matrix or vector
  single <- is.vector(po)
  
  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold <- !in_train
    
    if (single) {
      # Single vector case (oracle, semi-oracle)
      forest <- regression_forest(X[in_train, ], po[in_train])
    } else {
      # Matrix case (RF)
      forest <- regression_forest(X[in_train, ], po[in_train, fold])
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

# crossing fitting of nuisance parameters using RF
nuisance_rf <- function(X, Y, W, fold_indices, fold_pairs) {
  
  # double crossfitting across fold pairs
  cross_fits <- future_map(seq_along(fold_pairs), function(i) {
    fold_pair <- fold_pairs[[i]]
    in_train <- !(fold_indices %in% fold_pair)
    in_test <- !in_train
    
    Y.hat.model <- regression_forest(cbind(W[in_train], X[in_train, ]), Y[in_train])
    Y.hat.cf.model <- regression_forest(X[in_train, ], Y[in_train])
    W.hat.model <- regression_forest(X[in_train, ], W[in_train])
    
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
