###################
# functions for running DR-learner base CATE models - binary
# 
###################

require(grf, GenericML, dplyr, furrr, SuperLearner)

# Main function that runs all DR methods on a single dataset
run_all_dr_methods_optimized <- function(data, n_folds = 10, B = 200, workers = 5, 
                                         scenario = NULL, n = NULL, sl_lib = NULL) {
  
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
  
  # 1. DR Oracle ----
  if (!is.null(scenario) && !is.null(n)) {
    cat("Running DR Oracle...\n")
    results$oracle <- run_dr_oracle_optimized(data, X, Y, W, fmla_info, n_obs, fold_indices, 
                                              fold_list, scenario, n, n_folds, B, workers)
  }
  
  # 2. DR Semi-Oracle ----
  cat("Running DR Semi-Oracle...\n")
  results$semi_oracle <- run_dr_semi_oracle_optimized(data, X, Y, W, n_obs, fold_indices, 
                                                      fold_list, n_folds, B, workers)
  
  # 3. DR Random Forest ----
  cat("Running DR Random Forest...\n")
  results$random_forest <- run_dr_rf_optimized(data, X, Y, W, n_obs, fold_indices, 
                                               fold_list, fold_pairs, n_folds, B, workers)
  
  # 4. DR SuperLearner ----
  if (!is.null(sl_lib)) {
    cat("Running DR SuperLearner...\n")
    results$superlearner <- run_dr_superlearner_optimized(data, X, Y, W, n_obs, fold_indices, 
                                                          fold_list, fold_pairs, n_folds, 
                                                          B, workers, sl_lib)
  }
  
  return(results)
}

# Oracle-specific method (optimized)
run_dr_oracle_optimized <- function(data, X, Y, W, flma_info, n_obs, fold_indices, fold_list, 
                                    scenario, n, n_folds, B, workers) {
  
  # Load oracle formula info
  fmla <- parse(text = fmla_info[[1]])
  
  # Assign parameters from formula info
  list2env(fmla_info$params, envir = environment())
  
  # Oracle nuisance functions
  X_df <- as.data.frame(X)
  
  # Calculate oracle predictions
  W_temp <- rep(1, n_obs)
  Y1.hat <- plogis(eval(fmla, envir = list2env(c(list(W = W_temp), X_df))))
  
  W_temp <- rep(0, n_obs)
  Y0.hat <- plogis(eval(fmla, envir = list2env(c(list(W = W_temp), X_df))))
  
  Y.hat <- plogis(eval(fmla, envir = list2env(c(list(W = W), X_df))))
  W.hat <- rep(0.5, n_obs)
  
  # Calculate pseudo-outcomes
  cate <- Y1.hat - Y0.hat
  po <- cate + ((Y - Y.hat) * (W - W.hat)) / (W.hat * (1 - W.hat))
  
  # Estimate CATE and run analysis
  return(estimate_cate_and_analyze_optimized(X, Y, W, po, Y0.hat, W.hat, fold_indices, 
                                             fold_list, n_folds, B, workers))
}

# Semi-oracle method (optimized)
run_dr_semi_oracle_optimized <- function(data, X, Y, W, n_obs, fold_indices, fold_list, 
                                         n_folds, B, workers) {
  
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
  
  return(estimate_cate_and_analyze_optimized(X, Y, W, po, Y0.hat, W.hat, fold_indices, 
                                             fold_list, n_folds, B, workers))
}

# Random Forest method (optimized)
run_dr_rf_optimized <- function(data, X, Y, W, n_obs, fold_indices, fold_list, fold_pairs, 
                                n_folds, B, workers) {
  
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
    
    list(po = po, Y0.hat = Y0.hat, W.hat = W.hat, fold_pair = fold_pair, pair_index = i)
  }, .options = furrr_options(seed = TRUE))
  
  # Collate matrices using optimized approach
  po_matrix <- collate_predictions_optimized(fold_list, fold_pairs, fold_indices, cross_fits, "po")
  Y0.hat_matrix <- collate_predictions_optimized(fold_list, fold_pairs, fold_indices, cross_fits, "Y0.hat")
  W.hat_matrix <- collate_predictions_optimized(fold_list, fold_pairs, fold_indices, cross_fits, "W.hat")
  
  # Average predictions across folds
  po <- rowMeans(po_matrix, na.rm = TRUE)
  Y0.hat <- rowMeans(Y0.hat_matrix, na.rm = TRUE)
  W.hat <- rowMeans(W.hat_matrix, na.rm = TRUE)
  
  return(estimate_cate_and_analyze_optimized(X, Y, W, po, Y0.hat, W.hat, fold_indices, 
                                             fold_list, n_folds, B, workers, po_matrix))
}

# SuperLearner method (optimized)
run_dr_superlearner_optimized <- function(data, X, Y, W, n_obs, fold_indices, fold_list, 
                                          fold_pairs, n_folds, B, workers, sl_lib) {
  
  # Parallel cross-fitting with fold pairs
  cross_fits <- future_map(seq_along(fold_pairs), function(i) {
    fold_pair <- fold_pairs[[i]]
    in_train <- !(fold_indices %in% fold_pair)
    in_test <- which(!in_train)
    
    X_train <- X[in_train, ]
    X_W_train <- cbind(W = W[in_train], X_train)
    
    Y.hat.model <- SuperLearner(Y = Y[in_train], X = X_W_train, family = binomial(), 
                                SL.library = sl_lib, method = method.NNloglik())
    W.hat.model <- SuperLearner(W[in_train], X[in_train, ], family = binomial(), 
                                SL.library = sl_lib, method = method.NNloglik())
    
    X_test <- X[in_test, ]
    Y0.hat <- predict(Y.hat.model, newdata = cbind(W = 0, X_test))$pred
    Y1.hat <- predict(Y.hat.model, newdata = cbind(W = 1, X_test))$pred
    W.hat <- predict(W.hat.model, newdata = X_test)$pred
    
    # Handle failed W.hat model
    if (all(W.hat == 0)) {
      warning("All SuperLearner weights are zero for W.hat. Using mean(W).")
      W.hat <- rep(mean(W[in_test]), length(in_test))
    }
    
    # Trim extreme propensity scores
    W.hat[W.hat < 0.05] <- 0.05
    W.hat[W.hat > 0.95] <- 0.95
    
    W_test <- W[in_test]
    Y.hat <- W_test * Y1.hat + (1 - W_test) * Y0.hat
    
    cate <- Y1.hat - Y0.hat
    po <- cate + ((Y[in_test] - Y.hat) * (W_test - W.hat)) / (W.hat * (1 - W.hat))
    
    list(po = po, Y0.hat = Y0.hat, W.hat = W.hat, fold_pair = fold_pair, pair_index = i)
  }, .options = furrr_options(seed = TRUE))
  
  # Collate matrices
  po_matrix <- collate_predictions_optimized(fold_list, fold_pairs, fold_indices, cross_fits, "po")
  Y0.hat_matrix <- collate_predictions_optimized(fold_list, fold_pairs, fold_indices, cross_fits, "Y0.hat")
  W.hat_matrix <- collate_predictions_optimized(fold_list, fold_pairs, fold_indices, cross_fits, "W.hat")
  
  # Average predictions
  po <- rowMeans(po_matrix, na.rm = TRUE)
  Y0.hat <- rowMeans(Y0.hat_matrix, na.rm = TRUE)
  W.hat <- rowMeans(W.hat_matrix, na.rm = TRUE)
  
  return(estimate_cate_and_analyze_sl_optimized(X, Y, W, po, Y0.hat, W.hat, fold_indices, 
                                                fold_list, n_folds, B, workers, po_matrix, sl_lib))
}

# Optimized collate predictions function
collate_predictions_optimized <- function(fold_list, fold_pairs, fold_indices, cross_fits, target) {
  n_obs <- length(fold_indices)
  n_folds <- length(fold_list)
  
  # Initialize matrix
  result_matrix <- matrix(NA, nrow = n_obs, ncol = n_folds)
  
  # Fill matrix in parallel-friendly way
  for (i in seq_along(cross_fits)) {
    fold_pair <- cross_fits[[i]]$fold_pair
    predictions <- cross_fits[[i]][[target]]
    
    # Find which observations belong to this fold pair
    in_pair <- fold_indices %in% fold_pair
    
    # For each fold not in the current pair, assign predictions
    for (fold in fold_list) {
      if (!fold %in% fold_pair) {
        result_matrix[in_pair, fold] <- predictions
      }
    }
  }
  
  return(result_matrix)
}

# Optimized CATE estimation and analysis - RF based
estimate_cate_and_analyze_optimized <- function(X, Y, W, po, Y0.hat, W.hat, fold_indices, 
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
  BLP_tests <- future_map(seq_len(n_folds), function(fold) {
    in_fold <- fold_indices == fold
    tryCatch({
      BLP(Y[in_fold], W[in_fold], W.hat[in_fold], Y0.hat[in_fold], tau[in_fold])$coefficients[, c(1, 4)]
    }, error = function(e) {
      matrix(1, nrow = 4, ncol = 2)
    })
  }, .options = furrr_options(seed = TRUE))
  
  BLP_whole <- BLP(Y, W, W.hat, Y0.hat, tau)$coefficients[, c(1, 4)]
  
  # Get HTE p-values
  HTE_pval <- sapply(BLP_tests, function(x) {
    p_val <- x[4, 2]
    ifelse(is.na(p_val), 1, p_val)
  })
  
  # Confidence intervals
  ci_results <- get_confidence_intervals_optimized(X, Y, W, po, tau, fold_indices, fold_list, 
                                                   n_folds, B, workers, po_matrix)
  
  # TE-VIMs
  te_vims <- get_te_vims_optimized(X, po, tau, fold_indices, fold_list, n_folds, HTE_pval, 
                                   BLP_whole, po_matrix, workers)
  
  return(list(
    tau = ci_results,
    BLP_tests = BLP_tests,
    BLP_whole = BLP_whole,
    te_vims = te_vims,
    draws = ci_results$draws
  ))
}

# SuperLearner-specific analysis
estimate_cate_and_analyze_sl_optimized <- function(X, Y, W, po, Y0.hat, W.hat, fold_indices, 
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
  BLP_tests <- future_map(seq_len(n_folds), function(fold) {
    in_fold <- fold_indices == fold
    tryCatch({
      BLP(Y[in_fold], W[in_fold], W.hat[in_fold], Y0.hat[in_fold], tau[in_fold])$coefficients[, c(1, 4)]
    }, error = function(e) {
      matrix(1, nrow = 4, ncol = 2)
    })
  }, .options = furrr_options(seed = TRUE))
  
  BLP_whole <- BLP(Y, W, W.hat, Y0.hat, tau)$coefficients[, c(1, 4)]
  
  # Get HTE p-values
  HTE_pval <- sapply(BLP_tests, function(x) {
    p_val <- x[4, 2]
    ifelse(is.na(p_val), 1, p_val)
  })
  
  # Confidence intervals (SuperLearner version)
  ci_results <- get_confidence_intervals_sl_optimized(X, Y, W, po, tau, fold_indices, fold_list, 
                                                      n_folds, B, workers, po_matrix, sl_lib)
  
  # TE-VIMs (SuperLearner version)
  te_vims <- get_te_vims_sl_optimized(X, po, tau, fold_indices, fold_list, n_folds, HTE_pval, 
                                      BLP_whole, po_matrix, sl_lib, workers)
  
  return(list(
    tau = ci_results,
    BLP_tests = BLP_tests,
    BLP_whole = BLP_whole,
    te_vims = te_vims,
    draws = ci_results$draws
  ))
}

# Optimized confidence intervals 
get_confidence_intervals_optimized <- function(X, Y, W, po, tau, fold_indices, fold_list, 
                                               n_folds, B, workers, po_matrix = NULL) {
  
  # Bootstrap draws using future_map
  draws <- future_map(seq_len(B), function(b) {
    # Half samples
    half_samples <- lapply(fold_list, function(fold) {
      full <- sum(fold_indices == fold)
      half <- rep(FALSE, full)
      half[sample(1:full, floor(full/2), replace = FALSE)] <- TRUE
      half
    }) %>% unlist()
    
    # Compute half-sample CATEs
    tau_half <- lapply(fold_list, function(fold) {
      in_train <- half_samples & (fold_indices != fold)
      in_fold <- fold_indices == fold
      
      if (is.null(po_matrix)) {
        DR_rf <- regression_forest(X[in_train, ], po[in_train])
      } else {
        DR_rf <- regression_forest(X[in_train, ], po_matrix[in_train, fold])
      }
      
      predict(DR_rf, newdata = X[in_fold, ])$predictions
    }) %>% unlist() %>% unname()
    
    tau - tau_half
  }, .options = furrr_options(seed = TRUE))
  
  draws <- do.call(cbind, draws)
  
  # Confidence intervals
  lambda_hat <- apply(draws, 1, var)
  normalized <- abs(draws) / sqrt(lambda_hat)
  col_max <- apply(normalized, 2, max)
  S_star <- quantile(col_max, 0.975)
  
  data.frame(
    tau = tau,
    lb = tau - sqrt(lambda_hat) * S_star,
    ub = tau + sqrt(lambda_hat) * S_star,
    draws = I(draws)
  )
}

# Similarly optimized for SuperLearner
get_confidence_intervals_sl_optimized <- function(X, Y, W, po, tau, fold_indices, fold_list, 
                                                  n_folds, B, workers, po_matrix, sl_lib) {
  
  draws <- future_map(seq_len(B), function(b) {
    # Half samples
    half_samples <- lapply(fold_list, function(fold) {
      full <- sum(fold_indices == fold)
      half <- rep(FALSE, full)
      half[sample(1:full, floor(full/2), replace = FALSE)] <- TRUE
      half
    }) %>% unlist()
    
    # Compute half-sample CATEs
    tau_half <- lapply(fold_list, function(fold) {
      in_train <- half_samples & (fold_indices != fold)
      in_fold <- fold_indices == fold
      
      DR_rf <- SuperLearner(po_matrix[in_train, fold], X[in_train, ], 
                            family = gaussian(), SL.library = sl_lib)
      predict(DR_rf, newdata = X[in_fold, ])$pred
    }) %>% unlist() %>% unname()
    
    tau - tau_half
  }, .options = furrr_options(seed = TRUE))
  
  draws <- do.call(cbind, draws)
  
  # Confidence intervals
  lambda_hat <- apply(draws, 1, var)
  normalized <- abs(draws) / sqrt(lambda_hat)
  col_max <- apply(normalized, 2, max)
  S_star <- quantile(col_max, 0.975)
  
  data.frame(
    tau = tau,
    lb = tau - sqrt(lambda_hat) * S_star,
    ub = tau + sqrt(lambda_hat) * S_star,
    draws = I(draws)
  )
}

# Optimized TE-VIMs
get_te_vims_optimized <- function(X, po, tau, fold_indices, fold_list, n_folds, HTE_pval, 
                                  BLP_whole, po_matrix = NULL, workers) {
  
  if (!any(HTE_pval < 0.1) && BLP_whole[4, 2] >= 0.1) {
    return(NULL)
  }
  
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

# Optimized TE-VIMs for SuperLearner
get_te_vims_sl_optimized <- function(X, po, tau, fold_indices, fold_list, n_folds, HTE_pval, 
                                     BLP_whole, po_matrix, sl_lib, workers) {
  
  if (!any(HTE_pval < 0.1, na.rm = TRUE) && BLP_whole[4, 2] >= 0.1) {
    return(NULL)
  }
  
  n_obs <- nrow(X)
  X_df <- as.data.frame(X)
  covariates <- colnames(X_df)
  
  # Parallel computation of sub-taus for each covariate
  sub_taus_list <- future_map(seq_along(covariates), function(i) {
    new_X <- X_df %>% select(-!!covariates[i])
    
    # Adjust SL library for low-dimensional case
    sl_lib_adj <- if (length(covariates) == 2) {
      sl_lib[!sl_lib %in% "SL.glmnet"]
    } else {
      sl_lib
    }
    
    covariate_sub_taus <- rep(NA, n_obs)
    
    for (fold in seq_along(fold_list)) {
      in_train <- fold_indices != fold
      in_fold <- !in_train
      
      DR_sub <- SuperLearner(po_matrix[in_train, fold], new_X[in_train, ], 
                             SL.library = sl_lib_adj)
      covariate_sub_taus[in_fold] <- predict(DR_sub, newdata = new_X[in_fold, ])$pred
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