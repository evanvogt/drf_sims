###################
# Required packages
###################
require(coin)
require(grf)
require(SuperLearner)
require(future)
require(furrr)
require(dplyr)

###################
# functions for running all the CATE models - continuous
###################

# Main function that runs all methods on a single dataset
run_all_cate_methods <- function(data, n_folds = 10, B = 200, workers = 5, sl_lib = NULL, fmla_info = NULL) {
  
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
  
  # Results container
  results <- list()
  
  # First compute nuisance functions using DR Random Forest approach (most general)
  cat("Computing nuisance functions...\n")
  nuisance_results <- compute_nuisance_functions(X, Y, W, fold_indices, fold_pairs, workers)
  
  # 1. Causal Forest (reusing nuisance functions) ----
  cat("Running Causal Forest...\n")
  results$causal_forest <- run_causal_forest(X, Y, W, nuisance_results, fold_indices, 
                                             fold_list, n_folds, workers)
  
  # 2. DR Oracle ----
  if (!is.null(fmla_info)) {
    cat("Running DR Oracle...\n")
    results$dr_oracle <- run_dr_oracle(X, Y, W, fmla_info, fold_indices, 
                                       fold_list, n_folds, B, workers)
  }
  
  # 3. DR Semi-Oracle ----
  cat("Running DR Semi-Oracle...\n")
  results$dr_semi_oracle <- run_dr_semi_oracle(X, Y, W, fold_indices, 
                                               fold_list, n_folds, B, workers)
  
  # 4. DR Random Forest ----
  cat("Running DR Random Forest...\n")
  results$dr_random_forest <- run_dr_random_forest(X, Y, W, nuisance_results, fold_indices, 
                                                   fold_list, n_folds, B, workers)
  
  # 5. DR SuperLearner ----
  if (!is.null(sl_lib)) {
    cat("Running DR SuperLearner...\n")
    results$dr_superlearner <- run_dr_superlearner(X, Y, W, fold_indices, fold_list, 
                                                   fold_pairs, n_folds, B, workers, sl_lib)
  }
  
  return(results)
}

# Compute nuisance functions once (to be reused by CF and DR-RF)
compute_nuisance_functions <- function(X, Y, W, fold_indices, fold_pairs, workers) {
  
  # Parallel cross-fitting with fold pairs
  cross_fits <- future_map(seq_along(fold_pairs), function(i) {
    fold_pair <- fold_pairs[[i]]
    in_train <- !(fold_indices %in% fold_pair)
    in_test <- which(!in_train)
    
    Y.hat.model <- regression_forest(cbind(W[in_train], X[in_train, ]), Y[in_train])
    W.hat.model <- regression_forest(X[in_train, ], W[in_train])
    
    X_test <- X[in_test, ]
    
    Y0.hat <- predict(Y.hat.model, newdata = cbind(W = 0, X_test))$predictions
    Y1.hat <- predict(Y.hat.model, newdata = cbind(W = 1, X_test))$predictions
    W.hat <- predict(W.hat.model, newdata = X_test)$predictions
    
    W_test <- W[in_test]
    Y.hat <- W_test * Y1.hat + (1 - W_test) * Y0.hat
    
    cate <- Y1.hat - Y0.hat
    po <- cate + ((Y[in_test] - Y.hat) * (W_test - W.hat)) / (W.hat * (1 - W.hat))
    
    list(po = po, Y.hat = Y.hat, Y0.hat = Y0.hat, W.hat = W.hat, fold_pair = fold_pair)
  }, .options = furrr_options(seed = TRUE))
  
  # Collate predictions into matrices
  fold_list <- unique(fold_indices)
  po_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "po")
  Y.hat_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "Y.hat")
  Y0.hat_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "Y0.hat")
  W.hat_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "W.hat")
  
  return(list(
    po_matrix = po_matrix,
    Y.hat_matrix = Y.hat_matrix,
    Y0.hat_matrix = Y0.hat_matrix,
    W.hat_matrix = W.hat_matrix,
    po = rowMeans(po_matrix, na.rm = TRUE),
    Y0.hat = rowMeans(Y0.hat_matrix, na.rm = TRUE),
    W.hat = rowMeans(W.hat_matrix, na.rm = TRUE)
  ))
}

# Causal Forest method (reusing nuisance functions)
run_causal_forest <- function(X, Y, W, nuisance_results, fold_indices, fold_list, n_folds, workers) {
  
  n_obs <- nrow(X)
  
  # Get tau predictions using causal forest with pre-computed nuisance functions
  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold <- !in_train
    
    # Use nuisance predictions that didn't use the current fold
    forest <- causal_forest(X[in_train,], Y[in_train], W[in_train], 
                            nuisance_results$Y.hat_matrix[in_train, fold], 
                            nuisance_results$W.hat_matrix[in_train, fold])
    
    pred <- predict(forest, newdata = X[in_fold,], estimate.variance = TRUE)
    
    list(fold = fold, tau = pred$predictions, variance = pred$variance.estimates)
  }, .options = furrr_options(seed = TRUE))
  
  # Reconstruct tau vector
  tau <- rep(NA, n_obs)
  
  for (result in tau_results) {
    fold <- result$fold
    in_fold <- fold_indices == fold
    tau[in_fold] <- result$tau
  }
  
  # BLP and independence tests on the whole dataset only
  BLP_whole <- BLP(Y, W, nuisance_results$W.hat, nuisance_results$Y0.hat, tau)$coefficients[, c(1, 4)]
  independence_whole <- run_independence_test_whole(X, tau)
  
  return(list(
    tau = tau,
    BLP_whole = BLP_whole,
    independence_whole = independence_whole
  ))
}

# DR Oracle method
run_dr_oracle <- function(X, Y, W, fmla_info, fold_indices, fold_list, n_folds, B, workers) {
  
  n_obs <- nrow(X)
  
  # Assign parameters from formula info
  list2env(fmla_info$params, envir = environment())
  
  # Oracle nuisance functions
  X <- as.data.frame(X)
  
  # Calculate oracle predictions
  W_temp <- rep(1, n_obs)
  fmla <- parse(text = fmla_info$fmla)
  Y1.hat <- eval(fmla, envir = list2env(c(list(W = W_temp), X)))
  
  W_temp <- rep(0, n_obs)
  Y0.hat <- eval(fmla, envir = list2env(c(list(W = W_temp), X)))
  
  Y.hat <- eval(fmla, envir = list2env(c(list(W = W), X)))
  W.hat <- rep(0.5, n_obs)
  
  X <- as.matrix(X)
  # Calculate pseudo-outcomes
  cate <- Y1.hat - Y0.hat
  po <- cate + ((Y - Y.hat) * (W - W.hat)) / (W.hat * (1 - W.hat))
  
  # Estimate CATE and run analysis
  return(estimate_cate_and_analyze(X, Y, W, po, Y0.hat, W.hat, fold_indices, 
                                   fold_list, n_folds, B, workers))
}

# DR Semi-oracle method
run_dr_semi_oracle <- function(X, Y, W, fold_indices, fold_list, n_folds, B, workers) {
  
  n_obs <- nrow(X)
  W.hat <- rep(0.5, n_obs)
  
  # Parallel cross-fitting for Y model only
  cross_fits <- future_map(seq_len(n_folds), function(fold) {
    in_train <- !(fold_indices == fold)
    in_test <- which(!in_train)
    
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
  
  return(estimate_cate_and_analyze(X, Y, W, po, Y0.hat, W.hat, fold_indices, 
                                   fold_list, n_folds, B, workers))
}

# DR Random Forest method (reusing nuisance functions)
run_dr_random_forest <- function(X, Y, W, nuisance_results, fold_indices, fold_list, n_folds, B, workers) {
  
  return(estimate_cate_and_analyze(X, Y, W, nuisance_results$po, nuisance_results$Y0.hat, 
                                   nuisance_results$W.hat, fold_indices, fold_list, n_folds, 
                                   B, workers, nuisance_results$po_matrix))
}

# DR SuperLearner method
run_dr_superlearner <- function(X, Y, W, fold_indices, fold_list, fold_pairs, n_folds, B, workers, sl_lib) {
  X <- as.data.frame(X)
  on.exit({X <- as.matrix(X)})
  
  cross_fits <- future_map(seq_along(fold_pairs), function(i) {
    fold_pair <- fold_pairs[[i]]
    in_train <- !(fold_indices %in% fold_pair)
    in_test <- which(!in_train)
    
    X_train <- X[in_train, ]
    X_W_train <- cbind(W = W[in_train], X_train)
    
    # Fit outcome model
    Y.hat.model <- SuperLearner(Y = Y[in_train], X = X_W_train, SL.library = sl_lib)
    
    # Fit propensity model
    W.hat.model <- SuperLearner(W[in_train], X[in_train, ], family = binomial(), 
                                SL.library = sl_lib, method = method.NNloglik())
    
    X_test <- X[in_test, ]
    Y0.hat <- predict(Y.hat.model, newdata = cbind(W = 0, X_test))$pred
    Y1.hat <- predict(Y.hat.model, newdata = cbind(W = 1, X_test))$pred
    W.hat  <- predict(W.hat.model, newdata = X_test)$pred
    
    # --- NEW SAFEGUARDS ---
    # If SuperLearner fails to select learners and returns constant 0 values
    if (all(Y0.hat == 0) && all(Y1.hat == 0)) {
      warning("SuperLearner failed for Y.hat. Using mean(Y).")
      Y0.hat <- rep(mean(Y[in_train][W[in_train] == 0], na.rm = TRUE), length(in_test))
      Y1.hat <- rep(mean(Y[in_train][W[in_train] == 1], na.rm = TRUE), length(in_test))
    }
    
    if (all(W.hat == 0)) {
      warning("SuperLearner failed for W.hat. Using mean(W).")
      W.hat <- rep(mean(W[in_train], na.rm = TRUE), length(in_test))
    }
    
    # Trim extreme propensity scores
    W.hat[W.hat < 0.05] <- 0.05
    W.hat[W.hat > 0.95] <- 0.95
    
    # Pseudo-outcomes
    W_test <- W[in_test]
    Y.hat  <- W_test * Y1.hat + (1 - W_test) * Y0.hat
    cate   <- Y1.hat - Y0.hat
    po     <- cate + ((Y[in_test] - Y.hat) * (W_test - W.hat)) / (W.hat * (1 - W.hat))
    
    list(po = po, Y0.hat = Y0.hat, W.hat = W.hat, fold_pair = fold_pair)
  }, .options = furrr_options(seed = TRUE))
  
  # Collate matrices
  po_matrix    <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "po")
  Y0.hat_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "Y0.hat")
  W.hat_matrix  <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "W.hat")
  
  # Averages
  po    <- rowMeans(po_matrix, na.rm = TRUE)
  Y0.hat <- rowMeans(Y0.hat_matrix, na.rm = TRUE)
  W.hat  <- rowMeans(W.hat_matrix, na.rm = TRUE)
  
  return(estimate_cate_and_analyze_superlearner(
    X, Y, W, po, Y0.hat, W.hat, fold_indices, 
    fold_list, n_folds, B, workers, po_matrix, sl_lib
  ))
}
###################
# Helper Functions
###################

# Independence test on whole dataset only
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
  
  # Only tests on the whole dataset
  BLP_whole <- BLP(Y, W, W.hat, Y0.hat, tau)$coefficients[, c(1, 4)]
  independence_whole <- run_independence_test_whole(X, tau)
  
  return(list(
    tau = tau,
    BLP_whole = BLP_whole,
    independence_whole = independence_whole
  ))
}

# SuperLearner-specific analysis
estimate_cate_and_analyze_superlearner <- function(X, Y, W, po, Y0.hat, W.hat, fold_indices, 
                                                   fold_list, n_folds, B, workers, po_matrix, sl_lib) {
  
  n_obs <- nrow(X)
  
  # Parallel estimation of tau using SuperLearner
  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold <- !in_train
    
    forest <- SuperLearner(po_matrix[in_train, fold], X[in_train, ], 
                           family = gaussian(), SL.library = sl_lib)
    tau_pred <- predict(forest, newdata = X[in_fold, ])$pred
    # --------- SAFEGUARD -------
    if (all(tau_pred == 0)) {
      warning("SuperLearner returned all zero tau. Using mean(po_matrix) for fold ", fold, ".")
      tau_pred <- rep(mean(po_matrix[in_train, fold], na.rm = TRUE), sum(in_fold))
    }
    # ---------------------------
    list(fold = fold, predictions = tau_pred)
  }, .options = furrr_options(seed = TRUE))
  
  # Reconstruct tau vector
  tau <- rep(NA, n_obs)
  for (result in tau_results) {
    fold <- result$fold
    in_fold <- fold_indices == fold
    tau[in_fold] <- result$predictions
  }
  
  # Only tests on the whole dataset
  BLP_whole <- BLP(Y, W, W.hat, Y0.hat, tau)$coefficients[, c(1, 4)]
  independence_whole <- run_independence_test_whole(X, tau)
  
  return(list(
    tau = tau,
    BLP_whole = BLP_whole,
    independence_whole = independence_whole
  ))
}
