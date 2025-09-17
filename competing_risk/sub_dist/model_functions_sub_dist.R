###################
# Required packages
###################
require(coin)
require(grf)
require(future)
require(furrr)
require(dplyr)

###################
# Main function for running CATE models on competing risks survival data
###################
run_all_cate_methods_survival <- function(data, n_folds = 10, workers = 5, 
                                          horizon = NULL, fmla_info = NULL) {
  # Set up parallel processing plan
  metaplan <- plan(multisession, workers = workers)
  on.exit(plan(metaplan), add = TRUE)
  # Data preparation
  X <- as.matrix(data[, !(names(data) %in% c("Y", "D", "W"))])
  Y <- data$Y  # Event time
  D <- data$D  # Event type (0 = censored, 1 = event of interest, 2 = competing event)
  W <- data$W  # Treatment indicator
  n_obs <- nrow(X)
  # If no horizon specified, use maximum event time
  if (is.null(horizon)) {
    horizon <- max(Y[D != 0])  # Use max observed event time
    warning("No horizon specified. Using maximum observed event time: ", horizon)
  }
  # Create fold indices (common to all methods)
  fold_indices <- sort(seq(n_obs) %% n_folds) + 1
  fold_list <- unique(fold_indices)
  fold_pairs <- utils::combn(fold_list, 2, simplify = FALSE)
  # Results container
  results <- list()
  # 1. Causal Survival Forest ----
  cat("Running Causal Survival Forest...\n")
  results$causal_survival_forest <- run_causal_survival_forest(
    X, Y, D, W, fold_indices, fold_list, n_folds, horizon, workers
  )
  # 2. DR Oracle (if formula info provided) ----
  if (!is.null(fmla_info)) {
    cat("Running DR Oracle...\n")
    results$dr_oracle <- run_dr_oracle_survival(
      X, Y, D, W, fmla_info, fold_indices, fold_list, n_folds, horizon, workers
    )
  }
  # 3. DR Semi-Oracle ----
  cat("Running DR Semi-Oracle...\n")
  results$dr_semi_oracle <- run_dr_semi_oracle_survival(
    X, Y, D, W, fold_indices, fold_list, n_folds, horizon, workers
  )
  # 4. DR Survival Forest (using grf::survival_forest) ----
  cat("Running DR Survival Forest...\n")
  results$dr_survival_forest <- run_dr_survival_forest_grf(
    X, Y, D, W, fold_indices, fold_list, fold_pairs, n_folds, horizon, workers
  )
  return(results)
}

###################
# Helper function to handle censoring within a fold (Aalen-Johansen augmentation)
###################
handle_censoring_in_fold <- function(X, Y, D, W, horizon, in_train) {
  n_obs <- length(Y)
  # Create pseudo-outcomes for competing risks (Aalen-Johansen style augmentation)
  Y_pseudo <- Y
  Y_pseudo[D == 2] <- horizon + 1  # Move competing events beyond horizon
  D_pseudo <- as.integer(D == 1)   # Binary indicator for event of interest (1) or censored (0)
  # Check if there's censoring during follow-up in training set
  has_censoring <- any(D[in_train] == 0 & Y[in_train] < horizon)
  if (!has_censoring) {
    # No censoring during follow-up - simple case
    return(list(
      Y_pseudo = Y_pseudo,
      D_pseudo = D_pseudo,
      sample_weights = rep(1, n_obs),
      observed_events = rep(TRUE, n_obs),
      has_censoring = FALSE
    ))
  } else {
    # Handle censoring using IPCW - fit on training data only
    X_train <- X[in_train, , drop = FALSE]
    W_train <- W[in_train]
    Y_train <- Y[in_train]
    D_train <- D[in_train]
    # Fit time-to-censoring distribution using composite endpoints on training data
    sf_censor <- survival_forest(cbind(X_train, W_train), 
                                 Y = Y_train, 
                                 D = (D_train == 0))
    # Predict censoring probabilities for all observations
    sf_pred <- predict(sf_censor, 
                       newdata = cbind(X, W),
                       failure.times = sapply(Y, function(x) min(x, horizon)), 
                       prediction.times = "time")
    censoring_prob <- sf_pred$predictions
    # Identify observed events (not censored before horizon)
    observed_events <- (D != 0 | (D == 0 & Y >= horizon))
    # Calculate IPCW weights
    sample_weights <- rep(1, n_obs)
    sample_weights[observed_events] <- 1 / sapply(censoring_prob[observed_events], 
                                                  function(x) max(x, 1e-3))
    return(list(
      Y_pseudo = Y_pseudo,
      D_pseudo = D_pseudo,
      sample_weights = sample_weights,
      observed_events = observed_events,
      has_censoring = TRUE,
      censoring_prob = censoring_prob
    ))
  }
}

###################
# Method 1: Causal Survival Forest
###################
run_causal_survival_forest <- function(X, Y, D, W, fold_indices, 
                                       fold_list, n_folds, horizon, workers) {
  n_obs <- nrow(X)
  # Cross-validated tau estimation with fold-specific censoring weights
  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold <- !in_train
    # Handle censoring within this fold
    censoring_info <- handle_censoring_in_fold(X, Y, D, W, horizon, in_train)
    Y_pseudo <- censoring_info$Y_pseudo
    D_pseudo <- censoring_info$D_pseudo
    sample_weights <- censoring_info$sample_weights
    observed_events <- censoring_info$observed_events
    if (censoring_info$has_censoring) {
      # With censoring - use only observed events for training
      train_obs <- in_train & observed_events
      if (sum(train_obs) < 20) {
        return(list(fold = fold, tau = rep(NA, sum(in_fold))))
      }
      forest <- causal_survival_forest(
        X = X[train_obs, , drop = FALSE], 
        Y = Y_pseudo[train_obs], 
        D = D_pseudo[train_obs], 
        W = W[train_obs],
        sample.weights = sample_weights[train_obs],
        horizon = horizon
      )
    } else {
      # No censoring - use all training observations
      forest <- causal_survival_forest(
        X = X[in_train, , drop = FALSE], 
        Y = Y_pseudo[in_train], 
        D = D_pseudo[in_train], 
        W = W[in_train],
        horizon = horizon
      )
    }
    # Predict for test fold
    tau_fold <- predict(forest, newdata = X[in_fold, , drop = FALSE])$predictions
    list(fold = fold, tau = tau_fold)
  }, .options = furrr_options(seed = TRUE))
  # Reconstruct tau vector
  tau <- rep(NA, n_obs)
  for (result in tau_results) {
    fold <- result$fold
    in_fold <- fold_indices == fold
    tau[in_fold] <- result$tau
  }
  # For BLP tests, we need consistent pseudo outcomes and weights
  Y_pseudo <- Y
  Y_pseudo[D == 2] <- horizon + 1
  D_pseudo <- as.integer(D == 1)
  W_hat <- rep(0.5, n_obs)
  Y0_hat <- rep(mean(Y_pseudo[W == 0]), n_obs)
  BLP_tests <- run_blp_tests_survival(Y_pseudo, W, W_hat, Y0_hat, tau, 
                                      fold_indices, n_folds, workers)
  BLP_whole <- BLP_survival(Y_pseudo, W, W_hat, Y0_hat, tau)$coefficients[, c(1, 4)]
  independence_tests <- run_independence_tests_survival(X, tau, fold_indices, n_folds, workers)
  independence_whole <- run_independence_test_whole_survival(X, tau)
  return(list(
    tau = data.frame(tau = tau),
    BLP_tests = BLP_tests,
    BLP_whole = BLP_whole,
    independence_tests = independence_tests,
    independence_whole = independence_whole,
    method = "causal_survival_forest"
  ))
}

###################
# Method 2: DR Oracle for Survival (Weibull-based)
###################
run_dr_oracle_survival <- function(X, Y, D, W, fmla_info, fold_indices, 
                                   fold_list, n_folds, horizon, workers) {
  n_obs <- nrow(X)
  # Extract oracle parameters (based on your data generation)
  event1_params <- fmla_info$event1_params
  event2_params <- fmla_info$event2_params
  X_df <- as.data.frame(X)
  # Calculate true RMST for each individual under both treatments
  rmst_control <- numeric(n_obs)
  rmst_treat <- numeric(n_obs)
  for(i in 1:n_obs) {
    # Calculate scale parameters (prognostic effects)
    log_scale1 <- log(event1_params$scale_base) + 
      event1_params$b1_scale * X_df$X1[i] + 
      event1_params$b2_scale * X_df$X2[i]
    log_scale2 <- log(event2_params$scale_base) + 
      event2_params$b1_scale * X_df$X1[i] + 
      event2_params$b2_scale * X_df$X2[i]
    scale1 <- exp(log_scale1)
    scale2 <- exp(log_scale2)
    shape1_control <- event1_params$shape_base
    shape2_control <- event2_params$shape_base
    shape1_effect <- 0
    shape2_effect <- 0
    if (!is.null(event1_params$b3_shape) && "X3" %in% names(X_df)) {
      shape1_effect <- shape1_effect + event1_params$b3_shape * X_df$X3[i]
      shape2_effect <- shape2_effect + event2_params$b3_shape * X_df$X3[i]
    }
    if (!is.null(event1_params$b4_shape) && "X4" %in% names(X_df)) {
      shape1_effect <- shape1_effect + event1_params$b4_shape * X_df$X4[i]
      shape2_effect <- shape2_effect + event2_params$b4_shape * X_df$X4[i]
    }
    if (!is.null(event1_params$b34_shape) && all(c("X3", "X4") %in% names(X_df))) {
      shape1_effect <- shape1_effect + event1_params$b34_shape * X_df$X3[i] * X_df$X4[i]
      shape2_effect <- shape2_effect + event2_params$b34_shape * X_df$X3[i] * X_df$X4[i]
    }
    if (!is.null(event1_params$b45_shape) && all(c("X4", "X5") %in% names(X_df))) {
      shape1_effect <- shape1_effect + event1_params$b45_shape * X_df$X4[i] * X_df$X5[i]
      shape2_effect <- shape2_effect + event2_params$b45_shape * X_df$X4[i] * X_df$X5[i]
    }
    shape1_treat <- pmax(event1_params$shape_base + shape1_effect, 0.1)
    shape2_treat <- pmax(event2_params$shape_base + shape2_effect, 0.1)
    S1_control <- function(t) exp(-(t/scale1)^shape1_control)
    S2_control <- function(t) exp(-(t/scale2)^shape2_control)
    S_overall_control <- function(t) S1_control(t) * S2_control(t)
    rmst_control[i] <- integrate(S_overall_control, 0, horizon, 
                                 rel.tol = 1e-6)$value
    S1_treat <- function(t) exp(-(t/scale1)^shape1_treat)
    S2_treat <- function(t) exp(-(t/scale2)^shape2_treat)
    S_overall_treat <- function(t) S1_treat(t) * S2_treat(t)
    rmst_treat[i] <- integrate(S_overall_treat, 0, horizon,
                               rel.tol = 1e-6)$value
  }
  # Oracle pseudo-outcomes (true CATE)
  po <- rmst_treat - rmst_control
  # Estimate CATE using cross-validation
  tau <- estimate_cate_survival(X, po, fold_indices, fold_list, n_folds, workers)
  # For BLP tests
  Y_pseudo <- Y
  Y_pseudo[D == 2] <- horizon + 1
  W_hat <- rep(0.5, n_obs)
  Y0_hat <- rmst_control
  BLP_tests <- run_blp_tests_survival(Y_pseudo, W, W_hat, Y0_hat, tau, 
                                      fold_indices, n_folds, workers)
  BLP_whole <- BLP_survival(Y_pseudo, W, W_hat, Y0_hat, tau)$coefficients[, c(1, 4)]
  independence_tests <- run_independence_tests_survival(X, tau, fold_indices, n_folds, workers)
  independence_whole <- run_independence_test_whole_survival(X, tau)
  return(list(
    tau = data.frame(tau = tau),
    BLP_tests = BLP_tests,
    BLP_whole = BLP_whole,
    independence_tests = independence_tests,
    independence_whole = independence_whole,
    method = "dr_oracle_survival"
  ))
}

###################
# Method 3: DR Semi-Oracle for Survival
###################
run_dr_semi_oracle_survival <- function(X, Y, D, W, fold_indices, 
                                        fold_list, n_folds, horizon, workers) {
  n_obs <- nrow(X)
  W_hat <- rep(0.5, n_obs)
  # Cross-fitting for outcome model only
  cross_fits <- future_map(seq_len(n_folds), function(fold) {
    in_train <- !(fold_indices == fold)
    in_test <- which(!in_train)
    # Handle censoring within this fold
    censoring_info <- handle_censoring_in_fold(X, Y, D, W, horizon, in_train)
    Y_pseudo <- censoring_info$Y_pseudo
    D_pseudo <- censoring_info$D_pseudo
    sample_weights <- censoring_info$sample_weights
    observed_events <- censoring_info$observed_events
    if (censoring_info$has_censoring) {
      observed_train <- in_train & observed_events
      if (sum(observed_train) < 20) {
        return(list(po = rep(NA, length(in_test)), Y0_hat = rep(NA, length(in_test)), fold = fold))
      }
      Y_model <- survival_forest(
        cbind(W[observed_train], X[observed_train, , drop = FALSE]), 
        Y = Y_pseudo[observed_train], 
        D = D_pseudo[observed_train],
        sample.weights = sample_weights[observed_train]
      )
    } else {
      Y_model <- survival_forest(
        cbind(W[in_train], X[in_train, , drop = FALSE]), 
        Y = Y_pseudo[in_train], 
        D = D_pseudo[in_train]
      )
    }
    X_test <- X[in_test, , drop = FALSE]
    S0_hat <- predict(Y_model, newdata = cbind(W = 0, X_test),
                      failure.times = horizon)$predictions
    S1_hat <- predict(Y_model, newdata = cbind(W = 1, X_test),
                      failure.times = horizon)$predictions
    W_test <- W[in_test]
    S_hat <- W_test * S1_hat + (1 - W_test) * S0_hat
    cate <- S1_hat - S0_hat
    Y_test <- Y[in_test]
    D_test <- D[in_test]
    observed_survival <- (Y_test > horizon) | (D_test == 2)
    po <- cate + ((observed_survival - S_hat) * (W_test - 0.5)) / (0.5 * 0.5)
    list(po = po, Y0_hat = S0_hat, fold = fold)
  }, .options = furrr_options(seed = TRUE))
  po <- rep(NA, n_obs)
  Y0_hat <- rep(NA, n_obs)
  for (i in seq_along(cross_fits)) {
    fold <- cross_fits[[i]]$fold
    in_fold <- fold_indices == fold
    po[in_fold] <- cross_fits[[i]]$po
    Y0_hat[in_fold] <- cross_fits[[i]]$Y0_hat
  }
  tau <- estimate_cate_survival(X, po, fold_indices, fold_list, n_folds, workers)
  Y_pseudo <- Y
  Y_pseudo[D == 2] <- horizon + 1
  BLP_tests <- run_blp_tests_survival(Y_pseudo, W, W_hat, Y0_hat, tau, 
                                      fold_indices, n_folds, workers)
  BLP_whole <- BLP_survival(Y_pseudo, W, W_hat, Y0_hat, tau)$coefficients[, c(1, 4)]
  independence_tests <- run_independence_tests_survival(X, tau, fold_indices, n_folds, workers)
  independence_whole <- run_independence_test_whole_survival(X, tau)
  return(list(
    tau = data.frame(tau = tau),
    BLP_tests = BLP_tests,
    BLP_whole = BLP_whole,
    independence_tests = independence_tests,
    independence_whole = independence_whole,
    method = "dr_semi_oracle_survival"
  ))
}

###################
# Method 4: DR Survival Forest (using grf::survival_forest, not randomForestSRC)
###################
run_dr_survival_forest_grf <- function(X, Y, D, W, fold_indices, 
                                       fold_list, fold_pairs, n_folds, horizon, workers) {
  n_obs <- nrow(X)
  # Cross-fitting with fold pairs
  cross_fits <- future_map(seq_along(fold_pairs), function(i) {
    fold_pair <- fold_pairs[[i]]
    in_train <- !(fold_indices %in% fold_pair)
    in_test <- which(!in_train)
    # Handle censoring within this fold combination
    censoring_info <- handle_censoring_in_fold(X, Y, D, W, horizon, in_train)
    Y_pseudo <- censoring_info$Y_pseudo
    D_pseudo <- censoring_info$D_pseudo
    sample_weights <- censoring_info$sample_weights
    observed_events <- censoring_info$observed_events
    if (censoring_info$has_censoring) {
      observed_train <- in_train & observed_events
      if (sum(observed_train) < 20) {
        return(list(po = rep(NA, length(in_test)), 
                    Y0_hat = rep(NA, length(in_test)), 
                    W_hat = rep(0.5, length(in_test)),
                    fold_pair = fold_pair))
      }
      X_train <- cbind(W[observed_train], X[observed_train, , drop = FALSE])
      Y_pseudo_train <- Y_pseudo[observed_train]
      D_pseudo_train <- D_pseudo[observed_train]
      weights_train <- sample_weights[observed_train]
      Y_model <- survival_forest(
        X = as.matrix(X_train), 
        Y = Y_pseudo_train, 
        D = D_pseudo_train,
        sample.weights = weights_train
      )
    } else {
      X_train <- cbind(W[in_train], X[in_train, , drop = FALSE])
      Y_pseudo_train <- Y_pseudo[in_train]
      D_pseudo_train <- D_pseudo[in_train]
      Y_model <- survival_forest(
        X = as.matrix(X_train), 
        Y = Y_pseudo_train, 
        D = D_pseudo_train
      )
    }
    X_test_W0 <- cbind(0, X[in_test, , drop = FALSE])
    X_test_W1 <- cbind(1, X[in_test, , drop = FALSE])
    S0_hat <- predict(Y_model, newdata = as.matrix(X_test_W0), failure.times = horizon)$predictions
    S1_hat <- predict(Y_model, newdata = as.matrix(X_test_W1), failure.times = horizon)$predictions
    cate <- S1_hat - S0_hat
    W_test <- W[in_test]
    S_hat <- W_test * S1_hat + (1 - W_test) * S0_hat
    Y_test <- Y[in_test]
    D_test <- D[in_test]
    observed_survival <- (Y_test > horizon) | (D_test == 2)
    po <- cate + ((observed_survival - S_hat) * (W_test - 0.5)) / (0.25)
    list(po = po, Y0_hat = S0_hat, W_hat = rep(0.5, length(in_test)), fold_pair = fold_pair)
  }, .options = furrr_options(seed = TRUE))
  po_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "po")
  Y0_hat_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "Y0_hat")
  W_hat_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "W_hat")
  po <- rowMeans(po_matrix, na.rm = TRUE)
  Y0_hat <- rowMeans(Y0_hat_matrix, na.rm = TRUE)
  W_hat <- rowMeans(W_hat_matrix, na.rm = TRUE)
  tau <- estimate_cate_survival_with_matrix(X, po_matrix, fold_indices, fold_list, 
                                            n_folds, workers)
  Y_pseudo <- Y
  Y_pseudo[D == 2] <- horizon + 1
  BLP_tests <- run_blp_tests_survival(Y_pseudo, W, W_hat, Y0_hat, tau, 
                                      fold_indices, n_folds, workers)
  BLP_whole <- BLP_survival(Y_pseudo, W, W_hat, Y0_hat, tau)$coefficients[, c(1, 4)]
  independence_tests <- run_independence_tests_survival(X, tau, fold_indices, n_folds, workers)
  independence_whole <- run_independence_test_whole_survival(X, tau)
  return(list(
    tau = data.frame(tau = tau),
    BLP_tests = BLP_tests,
    BLP_whole = BLP_whole,
    independence_tests = independence_tests,
    independence_whole = independence_whole,
    method = "dr_survival_forest_grf"
  ))
}

###################
# Helper Functions
###################

estimate_cate_survival <- function(X, po, fold_indices, fold_list, n_folds, workers) {
  n_obs <- nrow(X)
  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold <- !in_train
    if (sum(!is.na(po[in_train])) < 10) {
      return(list(fold = fold, predictions = rep(NA, sum(in_fold))))
    }
    forest <- regression_forest(X[in_train, , drop = FALSE], po[in_train])
    tau_pred <- predict(forest, newdata = X[in_fold, , drop = FALSE])$predictions
    list(fold = fold, predictions = tau_pred)
  }, .options = furrr_options(seed = TRUE))
  tau <- rep(NA, n_obs)
  for (result in tau_results) {
    fold <- result$fold
    in_fold <- fold_indices == fold
    tau[in_fold] <- result$predictions
  }
  return(tau)
}

estimate_cate_survival_with_matrix <- function(X, po_matrix, fold_indices, fold_list, 
                                               n_folds, workers) {
  n_obs <- nrow(X)
  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold <- !in_train
    po_fold <- po_matrix[, fold]
    if (sum(!is.na(po_fold[in_train])) < 10) {
      return(list(fold = fold, predictions = rep(NA, sum(in_fold))))
    }
    forest <- regression_forest(X[in_train, , drop = FALSE], po_fold[in_train])
    tau_pred <- predict(forest, newdata = X[in_fold, , drop = FALSE])$predictions
    list(fold = fold, predictions = tau_pred)
  }, .options = furrr_options(seed = TRUE))
  tau <- rep(NA, n_obs)
  for (result in tau_results) {
    fold <- result$fold
    in_fold <- fold_indices == fold
    tau[in_fold] <- result$predictions
  }
  return(tau)
}

BLP_survival <- function(Y, W, W_hat, Y0_hat, tau) {
  n <- length(Y)
  complete_cases <- complete.cases(Y, W, W_hat, Y0_hat, tau)
  if (sum(complete_cases) < 10) {
    return(list(coefficients = matrix(NA, nrow = 4, ncol = 4)))
  }
  Y <- Y[complete_cases]
  W <- W[complete_cases]
  W_hat <- W_hat[complete_cases]
  Y0_hat <- Y0_hat[complete_cases]
  tau <- tau[complete_cases]
  tau_centered <- tau - mean(tau, na.rm = TRUE)
  X_blp <- cbind(1, W, W - W_hat, tau_centered * (W - W_hat))
  tryCatch({
    lm_fit <- lm(Y ~ X_blp - 1)
    return(summary(lm_fit))
  }, error = function(e) {
    return(list(coefficients = matrix(NA, nrow = 4, ncol = 4)))
  })
}

run_blp_tests_survival <- function(Y, W, W_hat, Y0_hat, tau, fold_indices, n_folds, workers) {
  future_map(seq_len(n_folds), function(fold) {
    in_fold <- fold_indices == fold
    tryCatch({
      BLP_survival(Y[in_fold], W[in_fold], W_hat[in_fold], Y0_hat[in_fold], tau[in_fold])$coefficients[, c(1, 4)]
    }, error = function(e) {
      matrix(1, nrow = 4, ncol = 2)
    })
  }, .options = furrr_options(seed = TRUE))
}

run_independence_tests_survival <- function(X, tau, fold_indices, n_folds, workers) {
  future_map(seq_len(n_folds), function(fold) {
    in_fold <- fold_indices == fold
    X_fold <- X[in_fold, ]
    tau_fold <- tau[in_fold]
    tryCatch({
      test_data <- data.frame(
        tau = tau_fold,
        X_fold
      )
      test_result <- coin::independence_test(
        tau ~ ., 
        data = test_data,
        distribution = "approximate",
        B = 9999
      )
      list(
        p_value = coin::pvalue(test_result),
        statistic = coin::statistic(test_result),
        method = "independence_test"
      )
    }, error = function(e) {
      list(
        p_value = 1,
        statistic = 0,
        method = "independence_test_failed"
      )
    })
  }, .options = furrr_options(seed = TRUE))
}

run_independence_test_whole_survival <- function(X, tau) {
  tryCatch({
    test_data <- data.frame(
      tau = tau,
      X
    )
    test_result <- coin::independence_test(
      tau ~ ., 
      data = test_data,
      distribution = "approximate",
      B = 9999
    )
    list(
      p_value = coin::pvalue(test_result),
      statistic = coin::statistic(test_result),
      method = "independence_test"
    )
  }, error = function(e) {
    list(
      p_value = 1,
      statistic = 0,
      method = "independence_test_failed"
    )
  })
}

collate_predictions <- function(fold_list, fold_pairs, fold_indices, cross_fits, target) {
  n_obs <- length(fold_indices)
  n_folds <- length(fold_list)
  result_matrix <- matrix(NA, nrow = n_obs, ncol = n_folds)
  for (i in seq_along(cross_fits)) {
    fold_pair <- cross_fits[[i]]$fold_pair
    predictions <- cross_fits[[i]][[target]]
    in_pair <- fold_indices %in% fold_pair
    for (fold in fold_list) {
      if (!fold %in% fold_pair) {
        result_matrix[in_pair, fold] <- predictions
        result_matrix[fold_indices == fold, fold] <- NA
      }
    }
  }
  return(result_matrix)
}
