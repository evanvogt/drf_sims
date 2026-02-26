##########
# title: CATE estimation functions - continuous outcome
##########
require(coin)
require(grf)
require(SuperLearner)
require(GenericML)
require(future)
require(furrr)
require(dplyr)

run_all_cate_methods <- function(data, n_folds = 10, sl_lib = NULL, fmla_info = NULL) {
  
  X <- as.matrix(data[, -c(1:2)]) 
  Y <- data$Y
  W <- data$W
  n_obs <- nrow(X)
  
  fold_indices <- sort(seq(n_obs) %% n_folds) + 1
  fold_list <- unique(fold_indices)
  fold_pairs <- utils::combn(fold_list, 2, simplify = FALSE)
  
  results <- list()
  
  # nuisance estimation (RF)
  cat("Computing nuisance functions...\n")
  nuisances_rf <- nuisance_rf(X, Y, W, fold_indices, fold_pairs)
  
  # Causal forest
  cat("Running Causal Forest...\n")
  results$causal_forest <- run_causal_forest(X, Y, W, nuisances_rf, fold_indices, fold_list)
  
  # DR random forest
  cat("Running DR Random Forest...\n")
  results$dr_random_forest <- run_dr_random_forest(X, Y, W, nuisances_rf, fold_indices, fold_list)
  
  # DR oracle
  if (!is.null(fmla_info)) {
    cat("Running DR Oracle...\n")
    results$dr_oracle <- run_dr_oracle(X, Y, W, fmla_info, fold_indices, fold_list)
  }
  
  # DR semi-oracle
  cat("Running DR Semi-Oracle...\n")
  results$dr_semi_oracle <- run_dr_semi_oracle(X, Y, W, fold_indices, fold_list)
  
  # DR SuperLearner
  if (!is.null(sl_lib)) {
    cat("Running DR SuperLearner...\n")
    X <- as.data.frame(X)
    nuisances_sl <- nuisance_sl(X, Y, W, fold_indices, fold_pairs, sl_lib)
    results$dr_superlearner <- run_dr_superlearner(X, Y, W, nuisances_sl, fold_indices, fold_list, fold_pairs, sl_lib)
    results$nuisances_sl <- nuisances_sl
  }
  
  results$nuisances_rf <- nuisances_rf
  results$fold_indices <- fold_indices
  
  return(results)
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
    po = rowMeans(po_matrix, na.rm = T),
    Y.hat_matrix = Y.hat_matrix,
    Y.hat = rowMeans(Y.hat_matrix, na.rm = T),
    Y.hat.cf_matrix = Y.hat.cf_matrix,
    Y.hat.cf = rowMeans(Y.hat.cf_matrix, na.rm = T),
    Y0.hat_matrix = Y0.hat_matrix,
    Y0.hat = rowMeans(Y0.hat_matrix, na.rm = T),
    W.hat_matrix = W.hat_matrix,
    W.hat = rowMeans(W.hat_matrix, na.rm = T)
  ))
}

# fit the causal forest
run_causal_forest <- function(X, Y, W, nuisances, fold_indices, fold_list) {
  n_obs <- nrow(X)
  
  # predict cates using nuisances
  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold <- !in_train
    
    forest <- causal_forest(X[in_train,], Y[in_train], W[in_train],
                            nuisances$Y.hat.cf_matrix[in_train, fold],
                            nuisances$W.hat_matrix[in_train, fold])
    
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
    BLP_whole = run_blp_whole(Y, W, nuisances$W.hat, nuisances$Y0.hat, tau),
    independence_cate = run_independence_test_whole(X, tau),
    independence_po = run_independence_test_whole(X, nuisances$po)
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

# DR learner with random forests
run_dr_random_forest <- function(X, Y, W, nuisances, fold_indices, fold_list) {
  tau <- stage_2_rf(X, nuisances$po_matrix, fold_indices, fold_list)
  
  return(list(
    tau = tau,
    BLP_whole = run_blp_whole(Y, W, nuisances$W.hat, nuisances$Y0.hat, tau),
    independence_cate = run_independence_test_whole(X, tau),
    independence_po = run_independence_test_whole(X, nuisances$po)
  ))
}

# DR RF oracle
run_dr_oracle <- function(X, Y, W, fmla_info, fold_indices, fold_list) {
  n_obs <- nrow(X)
  
  # Get parameters for oracle
  X <- as.data.frame(X)
  list2env(fmla_info$params, envir = environment())
  fmla <- parse(text = fmla_info$fmla)
  
  # calculate true values
  W_temp <- rep(1, n_obs)
  Y1.hat <- eval(fmla, envir = list2env(c(list(W = W_temp), X)))
  
  W_temp <- rep(0, n_obs)
  Y0.hat <- eval(fmla, envir = list2env(c(list(W = W_temp), X)))
  
  Y.hat <- eval(fmla, envir = list2env(c(list(W = W), X)))
  W.hat <- rep(0.5, n_obs)
  
  X <- as.matrix(X)
  
  # calculate pseudo outcomes
  cate <- Y1.hat - Y0.hat
  po <- cate + ((Y - Y.hat) * (W - W.hat)) / (W.hat * (1 - W.hat))
  
  tau <- stage_2_rf(X, po, fold_indices, fold_list)
  
  return(list(
    tau = tau,
    po = po,
    Y0.hat = Y0.hat,
    BLP_whole = run_blp_whole(Y, W, W.hat, Y0.hat, tau),
    independence_cate = run_independence_test_whole(X, tau),
    independence_po = run_independence_test_whole(X, po)
  ))
}

# DR RF semi-oracle
run_dr_semi_oracle <- function(X, Y, W, fold_indices, fold_list) {
  n_obs <- nrow(X)
  
  # set W.hat = 0.5
  W.hat <- rep(0.5, n_obs)
  
  # crossfit Y model only
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
    po <- cate + ((Y[in_test] - Y.hat) * (W_test - W.hat)) / (W.hat * (1 - W.hat))
    
    list(po = po, Y0.hat = Y0.hat, fold = fold)
  }, .options = furrr_options(seed = TRUE))
  
  # put results into a vector
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
    Y0.hat = Y0.hat,
    BLP_whole = run_blp_whole(Y, W, W.hat, Y0.hat, tau),
    independence_cate = run_independence_test_whole(X, tau),
    independence_po = run_independence_test_whole(X, po)
  ))
}

# crossfitting for nuisances with SL
nuisance_sl <- function(X, Y, W, fold_indices, fold_pairs, sl_lib) {
  
  # double crossfitting across fold pairs
  cross_fits <- future_map(seq_along(fold_pairs), function(i) {
    fold_pair <- fold_pairs[[i]]
    in_train <- !(fold_indices %in% fold_pair)
    in_test <- !in_train
    
    X_train <- X[in_train, ]
    X_W_train <- cbind(W = W[in_train], X_train)
    
    # train models + pretest SL libraries
    Y_lib <- pretest_superlearner(Y[in_train], X_W_train, sl_lib, gaussian())
    Y.hat.model <- SuperLearner(Y = Y[in_train], X = X_W_train, SL.library = Y_lib)
    
    W_lib <- pretest_superlearner(W[in_train], X_train, sl_lib, binomial())
    W.hat.model <- SuperLearner(W[in_train], X_train, family = binomial(), 
                                SL.library = W_lib, method = "method.NNloglik")
    
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
    po = rowMeans(po_matrix, na.rm = T),
    Y.hat_matrix = Y.hat_matrix,
    Y.hat = rowMeans(Y.hat_matrix, na.rm = T),
    Y0.hat_matrix = Y0.hat_matrix,
    Y0.hat = rowMeans(Y0.hat_matrix, na.rm = T),
    W.hat_matrix = W.hat_matrix,
    W.hat = rowMeans(W.hat_matrix, na.rm = T)
  ))
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

# second stage regression function (SuperLearner)
stage_2_sl <- function(X, po, fold_indices, fold_list, sl_lib) {
  n_obs <- nrow(X)
  
  # check if we have a po_matrix or vector
  single <- is.vector(po)
  
  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold <- !in_train
    
    if (single) {
      po_lib <- pretest_superlearner(po[in_train], X[in_train,], sl_lib, gaussian())
      po_model <- SuperLearner(po[in_train], X[in_train, ], family = gaussian(), SL.library = po_lib)
    } else {
      po_lib <- pretest_superlearner(po[in_train, fold], X[in_train,], sl_lib, gaussian())
      po_model <- SuperLearner(po[in_train, fold], X[in_train, ], family = gaussian(), SL.library = sl_lib)
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

run_dr_superlearner <- function(X, Y, W, nuisances, fold_indices, fold_list, fold_pairs, sl_lib) {

  tau <- stage_2_sl(X, nuisances$po_matrix, fold_indices, fold_list, sl_lib)
  
  return(list(
    tau = tau,
    BLP_whole = run_blp_whole(Y, W, nuisances$W.hat, nuisances$Y0.hat, tau),
    independence_cate = run_independence_test_whole(X, tau),
    independence_po = run_independence_test_whole(X, nuisances$po)
  ))
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
      teststat = "quadratic"
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