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
# Main CATE Methods Function
###################
run_all_cate_methods <- function(data, n_folds = 10, workers = 5, sl_lib = NULL, fmla_info = NULL) {
  metaplan <- plan(multisession, workers = workers)
  on.exit(plan(metaplan), add = TRUE)
  
  X <- as.matrix(data[, -c(1:2)]) 
  Y <- data$Y
  W <- data$W
  n_obs <- nrow(X)
  
  fold_indices <- sort(seq(n_obs) %% n_folds) + 1
  fold_list <- unique(fold_indices)
  fold_pairs <- utils::combn(fold_list, 2, simplify = FALSE)
  
  results <- list()
  
  # Compute nuisance functions (DR-RF) ----
  cat("Computing nuisance functions...\n")
  nuisance_results <- compute_nuisance_functions(X, Y, W, fold_indices, fold_pairs, workers)
  
  # 1. Causal Forest ----
  cat("Running Causal Forest...\n")
  results$causal_forest <- run_causal_forest(X, Y, W, nuisance_results, fold_indices, fold_list, workers)
  
  # 2. DR Oracle ----
  if (!is.null(fmla_info)) {
    cat("Running DR Oracle...\n")
    results$dr_oracle <- run_dr_oracle(X, Y, W, fmla_info, fold_indices, fold_list, workers)
  }
  
  # 3. DR Semi-Oracle ----
  cat("Running DR Semi-Oracle...\n")
  results$dr_semi_oracle <- run_dr_semi_oracle(X, Y, W, fold_indices, fold_list, workers)
  
  # 4. DR Random Forest ----
  cat("Running DR Random Forest...\n")
  results$dr_random_forest <- run_dr_random_forest(X, Y, W, nuisance_results, fold_indices, fold_list, workers)
  
  # 5. DR SuperLearner ----
  if (!is.null(sl_lib)) {
    cat("Running DR SuperLearner...\n")
    results$dr_superlearner <- run_dr_superlearner(X, Y, W, fold_indices, fold_list, fold_pairs, workers, sl_lib)
  }
  
  return(results)
}

###################
# Nuisance Functions
###################
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
    W.hat <- predict(W.hat.model, newdata = X_test)$predictions
    
    W_test <- W[in_test]
    Y.hat <- W_test * Y1.hat + (1 - W_test) * Y0.hat
    
    cate <- Y1.hat - Y0.hat
    po <- cate + ((Y[in_test] - Y.hat) * (W_test - W.hat)) / (W.hat * (1 - W.hat))
    
    list(po = po, Y0.hat = Y0.hat, W.hat = W.hat, fold_pair = fold_pair)
  }, .options = furrr_options(seed = TRUE))
  
  po_matrix    <- collate_predictions(unique(fold_indices), fold_pairs, fold_indices, cross_fits, "po")
  Y0.hat_matrix <- collate_predictions(unique(fold_indices), fold_pairs, fold_indices, cross_fits, "Y0.hat")
  W.hat_matrix  <- collate_predictions(unique(fold_indices), fold_pairs, fold_indices, cross_fits, "W.hat")
  
  return(list(
    po_matrix = po_matrix,
    Y0.hat_matrix = Y0.hat_matrix,
    W.hat_matrix = W.hat_matrix,
    po = rowMeans(po_matrix, na.rm = TRUE),
    Y0.hat = rowMeans(Y0.hat_matrix, na.rm = TRUE),
    W.hat = rowMeans(W.hat_matrix, na.rm = TRUE)
  ))
}

###################
# Estimation routines -- each returns tau and BLP and independence tests on full data only
###################
run_causal_forest <- function(X, Y, W, nuisance_results, fold_indices, fold_list, workers) {
  n_obs <- nrow(X)
  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold <- !in_train
    
    forest <- causal_forest(X[in_train,], Y[in_train], W[in_train],
                            nuisance_results$Y0.hat[in_train],
                            nuisance_results$W.hat[in_train])
    pred <- predict(forest, newdata = X[in_fold,])
    list(fold = fold, tau = pred$predictions)
  }, .options = furrr_options(seed = TRUE))
  
  tau <- rep(NA, n_obs)
  for (result in tau_results) {
    in_fold <- fold_indices == result$fold
    tau[in_fold] <- result$tau
  }
  
  return(list(
    tau = tau,
    BLP_whole = run_blp_whole(Y, W, nuisance_results$W.hat, nuisance_results$Y0.hat, tau),
    independence_whole = run_independence_test_whole(X, tau)
  ))
}

run_dr_oracle <- function(X, Y, W, fmla_info, fold_indices, fold_list, workers) {
  list2env(fmla_info$params, envir = environment())
  X <- as.data.frame(X)
  n_obs <- nrow(X)
  W_temp <- rep(1, n_obs)
  fmla <- parse(text = fmla_info$fmla)
  Y1.hat <- eval(fmla, envir = list2env(c(list(W = W_temp), X)))
  W_temp <- rep(0, n_obs)
  Y0.hat <- eval(fmla, envir = list2env(c(list(W = W_temp), X)))
  Y.hat <- eval(fmla, envir = list2env(c(list(W = W), X)))
  W.hat <- rep(0.5, n_obs)
  cate <- Y1.hat - Y0.hat
  po <- cate + ((Y - Y.hat) * (W - W.hat)) / (W.hat * (1 - W.hat))
  
  tau <- estimate_cate(X, po, fold_indices, fold_list, workers)
  return(list(
    tau = tau,
    BLP_whole = run_blp_whole(Y, W, W.hat, Y0.hat, tau),
    independence_whole = run_independence_test_whole(X, tau)
  ))
}

run_dr_semi_oracle <- function(X, Y, W, fold_indices, fold_list, workers) {
  n_obs <- nrow(X)
  W.hat <- rep(0.5, n_obs)
  
  cross_fits <- future_map(seq_len(length(fold_list)), function(fold) {
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
  
  po <- rep(NA, n_obs)
  Y0.hat <- rep(NA, n_obs)
  for (fit in cross_fits) {
    po[fold_indices == fit$fold] <- fit$po
    Y0.hat[fold_indices == fit$fold] <- fit$Y0.hat
  }
  
  tau <- estimate_cate(X, po, fold_indices, fold_list, workers)
  return(list(
    tau = tau,
    BLP_whole = run_blp_whole(Y, W, W.hat, Y0.hat, tau),
    independence_whole = run_independence_test_whole(X, tau)
  ))
}

run_dr_random_forest <- function(X, Y, W, nuisance_results, fold_indices, fold_list, workers) {
  tau <- estimate_cate(X, nuisance_results$po, fold_indices, fold_list, workers, nuisance_results$po_matrix)
  return(list(
    tau = tau,
    BLP_whole = run_blp_whole(Y, W, nuisance_results$W.hat, nuisance_results$Y0.hat, tau),
    independence_whole = run_independence_test_whole(X, tau)
  ))
}

run_dr_superlearner <- function(X, Y, W, fold_indices, fold_list, fold_pairs, workers, sl_lib) {
  X <- as.data.frame(X)
  
  cross_fits <- future_map(seq_along(fold_pairs), function(i) {
    fold_pair <- fold_pairs[[i]]
    in_train <- !(fold_indices %in% fold_pair)
    in_test <- which(!in_train)
    
    X_train <- X[in_train, ]
    X_W_train <- cbind(W = W[in_train], X_train)
    
    # pre-test base learners - failed algorithms will break the superlearner so must be removed manually to avoid error throwing
    Y_lib <- pretest_superlearner(Y[in_train], X_W_train, sl_lib, gaussian())
    Y.hat.model <- SuperLearner(Y[in_train], X_W_train, SL.library = Y_lib, method = "method.CC_LS")
    W_lib <- pretest_superlearner(W[in_train], X_train, sl_lib, binomial())
    W.hat.model <- SuperLearner(W[in_train], X_train, family = binomial(), SL.library = W_lib, method = "method.CC_LS")
    
    Y0.hat <- predict(Y.hat.model, newdata = cbind(W = 0, X[in_test,]))$pred
    Y1.hat <- predict(Y.hat.model, newdata = cbind(W = 1, X[in_test,]))$pred
    W.hat  <- predict(W.hat.model, newdata = X[in_test,])$pred
    
    
    W_test <- W[in_test]
    Y.hat <- W_test * Y1.hat + (1 - W_test) * Y0.hat
    cate <- Y1.hat - Y0.hat
    po <- cate + ((Y[in_test] - Y.hat) * (W_test - W.hat)) / (W.hat * (1 - W.hat))
    
    list(po = po, Y0.hat = Y0.hat, W.hat = W.hat, fold_pair = fold_pair)
  }, .options = furrr_options(seed = TRUE))
  
  po_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "po")
  Y0.hat_matrix <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "Y0.hat")
  W.hat_matrix  <- collate_predictions(fold_list, fold_pairs, fold_indices, cross_fits, "W.hat")
  
  po <- rowMeans(po_matrix, na.rm = TRUE)
  Y0.hat <- rowMeans(Y0.hat_matrix, na.rm = TRUE)
  W.hat  <- rowMeans(W.hat_matrix, na.rm = TRUE)
  
  W.hat[W.hat < 0.05] <- 0.05
  W.hat[W.hat > 0.95] <- 0.95
  
  tau <- c(rep(NA, n))
  for (fold in seq_along(fold_list)) {
    in_train <- fold_indices != fold
    in_fold <- !in_train
    tau_lib <- pretest_superlearner(po_matrix[in_train, fold], X[in_train,], sl_lib, gaussian())
    
    # If no learners work, assign mean
    if (length(tau_lib) == 0) {
      tau[in_fold] <- mean(po_matrix[in_train, fold], na.rm = TRUE)
      next
    }
    
    warning_flag <- FALSE
    sl_fit <- withCallingHandlers(
      {
        SuperLearner(po_matrix[in_train, fold], X[in_train,], family = gaussian(), SL.library = tau_lib)
      },
      warning = function(w) {
        if (any(grepl("All metalearner coefficients are zero, predictions will all be equal to 0", conditionMessage(w)))) {
          warning_flag <<- TRUE
          invokeRestart("muffleWarning")
        }
      }
    )
    # replace failed predictions with the mean
    if (warning_flag || is.null(sl_fit)) {
      tau[in_fold] <- mean(po_matrix[in_train, fold], na.rm = TRUE)
    } else {
      tau[in_fold] <- tryCatch(
        predict(sl_fit, newdata = X[in_fold, ])$pred,
        error = function(e) mean(po_matrix[in_train, fold], na.rm = TRUE)
      )
    }
  }
  
  
  return(list(
    tau = tau,
    BLP_whole = run_blp_whole(Y, W, W.hat, Y0.hat, tau),
    independence_whole = run_independence_test_whole(X, tau)
  ))
}

###################
# Helper Functions
###################
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

estimate_cate <- function(X, po, fold_indices, fold_list, workers, po_matrix = NULL, sl_lib = NULL) {
  n_obs <- nrow(X)
  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold <- !in_train
    
    if (is.null(sl_lib)) {
      if (is.null(po_matrix)) {
        forest <- regression_forest(X[in_train, ], po[in_train])
      } else {
        forest <- regression_forest(X[in_train, ], po_matrix[in_train, fold])
      }
      tau_pred <- predict(forest, newdata = X[in_fold, ])$predictions
    } else {
      forest <- SuperLearner(po_matrix[in_train, fold], X[in_train, ], family = gaussian(), SL.library = sl_lib)
      tau_pred <- predict(forest, newdata = X[in_fold, ])$pred
    }
    list(fold = fold, predictions = tau_pred)
  }, .options = furrr_options(seed = TRUE))
  
  tau <- rep(NA, n_obs)
  for (result in tau_results) {
    tau[fold_indices == result$fold] <- result$predictions
  }
  return(tau)
}

# Single BLP test on whole data
run_blp_whole <- function(Y, W, W.hat, Y0.hat, tau) {
  # Requires BLP function from your library
  BLP(Y, W, W.hat, Y0.hat, tau)$coefficients[, c(1, 4)]
}

# Single independence test on whole CATEs
run_independence_test_whole <- function(X, tau) {
  test_data <- data.frame(tau = tau, X)
  tryCatch({
    test_result <- coin::independence_test(
      tau ~ ., 
      data = test_data,
      teststat = "quadratic",
      distribution = "approximate"
    )
    list(
      p_value = coin::pvalue(test_result),
      statistic = coin::statistic(test_result),
      method = "independence_test"
    )
  }, error = function(e) {
    list(p_value = 1, statistic = 0, method = "independence_test_failed")
  })
}
