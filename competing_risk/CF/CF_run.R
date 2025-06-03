require(grf)
require(GenericML)
require(dplyr)
require(furrr)
require(future.apply)
CF_output <- function(data, n_folds, scenario) {
  # prepare the data ----
  X <- as.matrix(data[, -c(1:2)])  # Keeping only the covariates
  Y <- data$Y
  W <- data$W
  
  n <- dim(X)[1]
  
  # Create fold indices
  fold_indices <- sort(seq(n) %% n_folds) + 1
  fold_list <- unique(fold_indices)
  fold_pairs <- utils::combn(fold_list, 2, simplify = FALSE)  # Fold pairs
  
  # Y.hat, W.hat, po
  cross_fits <- lapply(fold_pairs, function(fold_pair) {
    in_train <- !(fold_indices %in% fold_pair) #train
    in_test <- which(!in_train)  #test
    
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
    
    list(po = po, Y.hat = Y.hat, Y0.hat = Y0.hat, W.hat = W.hat)
  })
  # matrix over predictions, excluding fold 1,2,3, etc
  po_matrix <- collate_predictions(seq_len(n_folds), fold_pairs, fold_indices, cross_fits, "po")
  Y.hat_matrix <- collate_predictions(seq_len(n_folds), fold_pairs, fold_indices, cross_fits, "Y.hat")
  Y0.hat_matrix <- collate_predictions(seq_len(n_folds), fold_pairs, fold_indices, cross_fits, "Y0.hat")
  W.hat_matrix <- collate_predictions(seq_len(n_folds), fold_pairs, fold_indices, cross_fits, "W.hat")
  
  # get tau predictions ----
  tau <- matrix(NA, n, 3)
  for (fold in seq_len(n_folds)) {
    in_train <- fold_indices != fold
    in_fold <- !in_train
    
    # take Y.hat and W.hat from the matrix, use predictions that didn't use the fold we want to predict
    forest <- causal_forest(X[in_train,], Y[in_train], W[in_train], 
                            Y.hat_matrix[in_train, fold], W.hat_matrix[in_train, fold]) 
    
    pred <- predict(forest, newdata = X[in_fold,], estimate.variance = TRUE)
    
    tau_est <- cbind(
      pred$predictions,
      pred$predictions - qnorm(0.975) * sqrt(pred$variance.estimates),
      pred$predictions + qnorm(0.975) * sqrt(pred$variance.estimates)
    )
    tau[in_fold, ] <- tau_est
  }
  tau <- as.data.frame(tau)
  colnames(tau) <- c("tau", "lb", "ub")
  
  # BLP tests ----
  Y.hat <- Y.hat_matrix %>% rowMeans(, na.rm = T)
  Y0.hat <- Y0.hat_matrix %>% rowMeans(, na.rm = T)
  W.hat <- W.hat_matrix %>% rowMeans(, na.rm = T)
  
  BLP_tests <- lapply(seq_len(n_folds), function(fold) {
    in_fold <- fold_indices == fold
    blp_test <- BLP(Y[in_fold], W[in_fold], W.hat[in_fold], Y0.hat[in_fold], tau$tau[in_fold])$coefficients[,c(1,4)]
    return(blp_test)
  })
  
  #BLP on the whole dataset
  BLP_whole <- BLP(Y, W, W.hat, Y0.hat, tau$tau)$coefficients[,c(1,4)]
  
  # make the te-vims conditional on there being a HTE signal
  HTE_pval <- lapply(BLP_tests, function(x) {
    p_val <- x[4,2]
  }) %>% unlist()
  
  print(paste0("BLP test result: ", min(HTE_pval, BLP_whole[4,2])))
  
  # TE-vims ----
  te_vims <- NULL
  if (any( HTE_pval < 0.1) | BLP_whole[4,2] < 0.1) {
    print("HTE signal - running TE_VIMS")
    # Rerunning the forests with leave-one-out method
    covariates <- colnames(X)
    sub_taus <- matrix(nrow = n, ncol = length(covariates))
    colnames(sub_taus) <- covariates
    
    compute_fold <- function(i) {
      cov <- covariates[i]
      new_X <- as.matrix(X[, -i])
      sub_taus_col <- numeric(n)
      
      for (fold in seq_along(fold_list)) {
        in_train <- fold_indices != fold
        in_fold <- !in_train
        new_forest <- causal_forest(as.matrix(new_X[in_train, ]), Y[in_train], W[in_train], Y.hat_matrix[in_train, fold], W.hat_matrix[in_train, fold])
        sub_taus_col <- predict(new_forest, newdata = as.matrix(new_X[in_fold, ]))$predictions
      }
      return(sub_taus_col)
    }
    
    results <- future_lapply(seq_along(covariates), compute_fold, future.seed = T)
    sub_taus <- do.call(cbind, results)
    colnames(sub_taus) <- covariates
    
    # Compute TE-VIMs
    po <- po_matrix %>% rowMeans(, na.rm = T)
    ate <- sum(po) / n
    
    r_ate <- (po - ate)^2
    r_tau <- (po - tau)^2
    
    te_vims <- apply(sub_taus, 2, function(sub_tau) {
      r_subtau <- (po - sub_tau)^2
      
      # evaluate TE-VIM (Theta_s in the paper)
      tevim <- sum(r_subtau - r_tau) / n
      
      infl <- r_subtau - r_tau - tevim
      std_err <- sqrt(sum(infl^2)) / n
      
      return(list(tevim = tevim, std_err = std_err))
    }) %>% simplify2array()
    
    te_vims <- as.data.frame(te_vims)
    
  }
  # output all results ---
  return(list(tau = tau, BLP_tests = BLP_tests, BLP_whole = BLP_whole, te_vims = te_vims))
}
