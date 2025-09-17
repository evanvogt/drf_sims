###################
# Required packages
###################
require(coin)
require(grf)
require(future)
require(furrr)
require(dplyr)
require(GenericML)

###################
# functions for running CATE models
###################

# Main function that runs Causal Forest and DR Random Forest
run_all_cate_methods <- function(data, n_folds = 10, workers = 5) {
  
  # Set up parallel processing plan
  metaplan <- plan(multisession, workers = workers)
  on.exit(plan(metaplan), add = TRUE)
  
  # Data preparation
  X <- as.matrix(data[, -c(1:2)]) 
  Y <- data$Y
  W <- data$W
  n_obs <- nrow(X)
  
  # Create fold indices (common to all methods)
  fold_indices <- sort(seq(n_obs) %% n_folds) + 1
  fold_list <- unique(fold_indices)
  fold_pairs <- utils::combn(fold_list, 2, simplify = FALSE)
  
  results <- list()
  
  # Compute nuisance functions
  cat("Computing nuisance functions...\n")
  nuisance_results <- compute_nuisance_functions(X, Y, W, fold_indices, fold_pairs, workers)
  
  # Causal Forest
  cat("Running Causal Forest...\n")
  results$causal_forest <- run_causal_forest(
    X, Y, W, nuisance_results, fold_indices, fold_list, n_folds, workers
  )
  
  # DR Random Forest
  cat("Running DR Random Forest...\n")
  results$dr_random_forest <- run_dr_random_forest(
    X, Y, W, nuisance_results, fold_indices, fold_list, n_folds, workers
  )
  
  return(results)
}

# Compute nuisance functions once (for CF & DR-RF)
compute_nuisance_functions <- function(X, Y, W, fold_indices, fold_pairs, workers) {
  
  cross_fits <- future_map(seq_along(fold_pairs), function(i) {
    fold_pair <- fold_pairs[[i]]
    in_train <- !(fold_indices %in% fold_pair)
    in_test <- which(!in_train)
    
    Y.hat.model <- regression_forest(cbind(W[in_train], X[in_train, ]), Y[in_train])
    W.hat.model <- regression_forest(X[in_train, ], W[in_train])
    
    X_test <- X[in_test, ]
    Y0.hat <- predict(Y.hat.model, newdata = cbind(W = 0, X_test))$predictions
    Y1.hat <- predict(Y.hat.model, newdata = cbind(W = 1, X_test))$predictions
    W.hat  <- predict(W.hat.model, newdata = X_test)$predictions
    
    W_test <- W[in_test]
    Y.hat  <- W_test * Y1.hat + (1 - W_test) * Y0.hat
    cate   <- Y1.hat - Y0.hat
    po     <- cate + ((Y[in_test] - Y.hat) * (W_test - W.hat)) / (W.hat * (1 - W.hat))
    
    list(po = po, Y.hat = Y.hat, Y0.hat = Y0.hat, W.hat = W.hat, fold_pair = fold_pair)
  }, .options = furrr_options(seed = TRUE))
  
  fold_list <- unique(fold_indices)
  po_matrix   <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "po")
  Y.hat_matrix<- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "Y.hat")
  Y0.hat_matrix<-collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "Y0.hat")
  W.hat_matrix<- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "W.hat")
  
  return(list(
    po_matrix   = po_matrix,
    Y.hat_matrix= Y.hat_matrix,
    Y0.hat_matrix=Y0.hat_matrix,
    W.hat_matrix=W.hat_matrix,
    po        = rowMeans(po_matrix, na.rm = TRUE),
    Y0.hat    = rowMeans(Y0.hat_matrix, na.rm = TRUE),
    W.hat     = rowMeans(W.hat_matrix, na.rm = TRUE)
  ))
}

# Causal Forest method with TE-VIMs
run_causal_forest <- function(X, Y, W, nuisance_results, fold_indices, fold_list, n_folds, workers) {
  n_obs <- nrow(X)
  
  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold <- !in_train
    
    forest <- causal_forest(X[in_train,], Y[in_train], W[in_train], 
                            nuisance_results$Y.hat_matrix[in_train, fold], 
                            nuisance_results$W.hat_matrix[in_train, fold])
    
    pred <- predict(forest, newdata = X[in_fold,], estimate.variance = FALSE)
    list(fold = fold, tau = pred$predictions)
  }, .options = furrr_options(seed = TRUE))
  
  tau <- rep(NA, n_obs)
  for (result in tau_results) {
    tau[fold_indices == result$fold] <- result$tau
  }
  
  BLP_tests <- run_blp_tests(Y, W, nuisance_results$W.hat, nuisance_results$Y0.hat, tau, fold_indices, n_folds, workers)
  BLP_whole <- BLP(Y, W, nuisance_results$W.hat, nuisance_results$Y0.hat, tau)$coefficients[, c(1, 4)]
  
  independence_tests <- run_independence_tests(X, tau, fold_indices, n_folds, workers)
  independence_whole <- run_independence_test_whole(X, tau)
  
  # Get HTE p-values
  HTE_pval <- sapply(BLP_tests, function(x) {
    p_val <- x[4, 2]
    ifelse(is.na(p_val), 1, p_val)
  })
  
  # TE-VIMs (Causal Forest version)
  te_vims <- get_te_vims_causal_forest(
    X, Y, W, nuisance_results, tau, fold_indices, fold_list, n_folds, HTE_pval, BLP_whole, workers
  )
  
  return(list(
    tau = tau,
    BLP_tests = BLP_tests,
    BLP_whole = BLP_whole,
    independence_tests = independence_tests,
    independence_whole = independence_whole,
    te_vims = te_vims
  ))
}

# DR Random Forest method with TE-VIMs
run_dr_random_forest <- function(X, Y, W, nuisance_results, fold_indices, fold_list, n_folds, workers) {
  n_obs <- nrow(X)
  
  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold <- !in_train
    
    forest <- regression_forest(X[in_train, ], nuisance_results$po_matrix[in_train, fold])
    tau_pred <- predict(forest, newdata = X[in_fold, ])$predictions
    
    list(fold = fold, predictions = tau_pred)
  }, .options = furrr_options(seed = TRUE))
  
  tau <- rep(NA, n_obs)
  for (result in tau_results) {
    tau[fold_indices == result$fold] <- result$predictions
  }
  
  BLP_tests <- run_blp_tests(Y, W, nuisance_results$W.hat, nuisance_results$Y0.hat, tau, fold_indices, n_folds, workers)
  BLP_whole <- BLP(Y, W, nuisance_results$W.hat, nuisance_results$Y0.hat, tau)$coefficients[, c(1, 4)]
  
  independence_tests <- run_independence_tests(X, tau, fold_indices, n_folds, workers)
  independence_whole <- run_independence_test_whole(X, tau)
  
  # Get HTE p-values
  HTE_pval <- sapply(BLP_tests, function(x) {
    p_val <- x[4, 2]
    ifelse(is.na(p_val), 1, p_val)
  })
  
  # TE-VIMs (DR Random Forest version)
  te_vims <- get_te_vims(
    X, nuisance_results$po, tau, fold_indices, fold_list, n_folds, HTE_pval, BLP_whole, 
    nuisance_results$po_matrix, workers
  )
  
  return(list(
    tau = tau,
    BLP_tests = BLP_tests,
    BLP_whole = BLP_whole,
    independence_tests = independence_tests,
    independence_whole = independence_whole,
    te_vims = te_vims
  ))
}

###################
# Helper Functions
###################

# BLP tests
run_blp_tests <- function(Y, W, W.hat, Y0.hat, tau, fold_indices, n_folds, workers) {
  future_map(seq_len(n_folds), function(fold) {
    in_fold <- fold_indices == fold
    tryCatch({
      BLP(Y[in_fold], W[in_fold], W.hat[in_fold], Y0.hat[in_fold], tau[in_fold])$coefficients[, c(1, 4)]
    }, error = function(e) {
      matrix(1, nrow = 4, ncol = 2)
    })
  }, .options = furrr_options(seed = TRUE))
}

# Independence tests
run_independence_tests <- function(X, tau, fold_indices, n_folds, workers) {
  future_map(seq_len(n_folds), function(fold) {
    in_fold <- fold_indices == fold
    X_fold <- X[in_fold, ]
    tau_fold <- tau[in_fold]
    
    tryCatch({
      # Convert to data frame for coin package
      test_data <- data.frame(
        tau = tau_fold,
        X_fold
      )
      
      # Run independence test using coin package
      # Tests if tau is independent of all covariates X
      test_result <- coin::independence_test(
        tau ~ ., 
        data = test_data,
        distribution = "approximate",
        B = 9999  # number of Monte Carlo replications
      )
      
      # Extract p-value and test statistic
      list(
        p_value = coin::pvalue(test_result),
        statistic = coin::statistic(test_result),
        method = "independence_test"
      )
      
    }, error = function(e) {
      # Return default values if test fails
      list(
        p_value = 1,
        statistic = 0,
        method = "independence_test_failed"
      )
    })
  }, .options = furrr_options(seed = TRUE))
}

# Independence test on whole dataset
run_independence_test_whole <- function(X, tau) {
  tryCatch({
    # Convert to data frame for coin package
    test_data <- data.frame(
      tau = tau,
      X
    )
    
    # Run independence test using coin package
    test_result <- coin::independence_test(
      tau ~ ., 
      data = test_data,
      distribution = "approximate",
      B = 9999  # number of Monte Carlo replications
    )
    
    # Extract results
    list(
      p_value = coin::pvalue(test_result),
      statistic = coin::statistic(test_result),
      method = "independence_test"
    )
    
  }, error = function(e) {
    # Return default values if test fails
    list(
      p_value = 1,
      statistic = 0,
      method = "independence_test_failed"
    )
  })
}

# Collate predictions function (uses existing function from your scripts)
collate_predictions <- function(fold_list, fold_pairs, fold_indices, cross_fits, target) {
  n_obs <- length(fold_indices)
  n_folds <- length(fold_list)
  
  # Initialize matrix
  result_matrix <- matrix(NA, nrow = n_obs, ncol = n_folds)
  
  # Fill matrix
  for (i in seq_along(cross_fits)) {
    fold_pair <- cross_fits[[i]]$fold_pair
    predictions <- cross_fits[[i]][[target]]
    
    # Find which observations belong to this fold pair
    in_pair <- fold_indices %in% fold_pair
    
    # For each fold not in the current pair, assign predictions
    for (fold in fold_list) {
      if (!fold %in% fold_pair) {
        result_matrix[in_pair, fold] <- predictions
        result_matrix[fold_indices == fold, fold] <- NA
      }
    }
  }
  
  return(result_matrix)
}

# CATE estimation and analysis - main version
estimate_cate_and_analyze <- function(X, Y, W, po, Y0.hat, W.hat, fold_indices, 
                                      fold_list, n_folds, B, workers, po_matrix = NULL) {
  
  n_obs <- nrow(X)
  
  # Parallel estimation of tau using cross-validation
  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold <- !in_train
    
    if (is.null(po_matrix)) {
      # Single vector case (oracle, semi-oracle)
      forest <- regression_forest(X[in_train, ], po[in_train])
    } else {
      # Matrix case (RF)
      forest <- regression_forest(X[in_train, ], po_matrix[in_train, fold])
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
  
  # BLP tests 
  BLP_tests <- run_blp_tests(Y, W, W.hat, Y0.hat, tau, fold_indices, n_folds, workers)
  BLP_whole <- BLP(Y, W, W.hat, Y0.hat, tau)$coefficients[, c(1, 4)]
  
  # Independence tests
  independence_tests <- run_independence_tests(X, tau, fold_indices, n_folds, workers)
  independence_whole <- run_independence_test_whole(X, tau)
  
  # Get HTE p-values
  HTE_pval <- sapply(BLP_tests, function(x) {
    p_val <- x[4, 2]
    ifelse(is.na(p_val), 1, p_val)
  })
  
  # Confidence intervals
  ci_results <- get_confidence_intervals(X, Y, W, po, tau, fold_indices, fold_list, 
                                         n_folds, B, workers, po_matrix)
  
  # TE-VIMs
  te_vims <- get_te_vims(X, po, tau, fold_indices, fold_list, n_folds, HTE_pval, 
                         BLP_whole, po_matrix, workers)
  
  return(list(
    tau = ci_results,
    BLP_tests = BLP_tests,
    BLP_whole = BLP_whole,
    independence_tests = independence_tests,
    independence_whole = independence_whole,
    te_vims = te_vims,
    draws = ci_results$draws
  ))
}


# TE-VIMs
get_te_vims <- function(X, po, tau, fold_indices, fold_list, n_folds, HTE_pval, 
                        BLP_whole, po_matrix = NULL, workers) {
  n_obs <- nrow(X)
  covariates <- colnames(X)
  
  # Parallel computation of sub-taus for each covariate
  sub_taus_list <- future_map(seq_along(covariates), function(i) {
    new_X <- as.matrix(X[, -i])
    covariate_sub_taus <- rep(NA, n_obs)
    
    for (fold in seq_along(fold_list)) {
      in_train <- fold_indices != fold
      in_fold <- !in_train
      
      if (is.null(po_matrix)) {
        DR_sub <- regression_forest(new_X[in_train, ], po[in_train])
      } else {
        DR_sub <- regression_forest(new_X[in_train, ], po_matrix[in_train, fold])
      }
      
      covariate_sub_taus[in_fold] <- predict(DR_sub, newdata = new_X[in_fold, ])$predictions
    }
    
    covariate_sub_taus
  }, .options = furrr_options(seed = TRUE))
  
  # Convert to matrix
  sub_taus <- do.call(cbind, sub_taus_list)
  colnames(sub_taus) <- covariates
  
  # Compute TE-VIMs
  ate <- mean(po)
  r_ate <- (po - ate)^2
  r_tau <- (po - tau)^2
  
  te_vims <- apply(sub_taus, 2, function(sub_tau) {
    r_subtau <- (po - sub_tau)^2
    tevim <- sum(r_subtau - r_tau) / n_obs
    infl <- r_subtau - r_tau - tevim
    std_err <- sqrt(sum(infl^2)) / n_obs
    list(tevim = tevim, std_err = std_err)
  }) %>% simplify2array()
  
  as.data.frame(te_vims)
}

# TE-VIMs for Causal Forest
get_te_vims_causal_forest <- function(X, Y, W, nuisance_results, tau, fold_indices, 
                                      fold_list, n_folds, HTE_pval, BLP_whole, workers) {
  n_obs <- nrow(X)
  covariates <- colnames(X)
  po <- nuisance_results$po
  
  # Parallel computation of sub-taus for each covariate
  sub_taus_list <- future_map(seq_along(covariates), function(i) {
    new_X <- as.matrix(X[, -i])
    covariate_sub_taus <- rep(NA, n_obs)
    
    for (fold in seq_along(fold_list)) {
      in_train <- fold_indices != fold
      in_fold <- !in_train
      
      # Use causal forest with reduced covariates
      forest <- causal_forest(new_X[in_train, ], Y[in_train], W[in_train], 
                              nuisance_results$Y.hat_matrix[in_train, fold], 
                              nuisance_results$W.hat_matrix[in_train, fold])
      
      covariate_sub_taus[in_fold] <- predict(forest, newdata = new_X[in_fold, ])$predictions
    }
    
    covariate_sub_taus
  }, .options = furrr_options(seed = TRUE))
  
  # Convert to matrix
  sub_taus <- do.call(cbind, sub_taus_list)
  colnames(sub_taus) <- covariates
  
  # Compute TE-VIMs
  ate <- mean(po)
  r_ate <- (po - ate)^2
  r_tau <- (po - tau)^2
  
  te_vims <- apply(sub_taus, 2, function(sub_tau) {
    r_subtau <- (po - sub_tau)^2
    tevim <- sum(r_subtau - r_tau) / n_obs
    infl <- r_subtau - r_tau - tevim
    std_err <- sqrt(sum(infl^2)) / n_obs
    list(tevim = tevim, std_err = std_err)
  }) %>% simplify2array()
  
  as.data.frame(te_vims)
}
