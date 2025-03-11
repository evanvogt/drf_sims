require(grf, GenericML, dplyr, furrr)
T_RF_output <- function(data, n_folds, scenario, B, workers) {
  X <- as.matrix(data[, -c(1:2)])  # Keeping only the covariates
  Y <- data$Y
  W <- data$W
  
  n <- dim(X)[1]
  
  # Create fold indices
  fold_indices <- sort(seq(n) %% n_folds) + 1
  
  # Perform cross-fits
  cross_fits <- lapply(seq_len(n_folds), function(fold) {
    in_train <- !(fold_indices == fold) #train
    in_test <- which(!in_train)  #test
    train0 <- in_train & W==0 # train control
    train1 <- in_train & W==1 # train treated
    
    Y0.hat.model <- regression_forest(X[train0, ], Y[train0])
    Y1.hat.model <- regression_forest(X[train1, ], Y[train1])
    W.hat.model <- regression_forest(X[in_train, ], W[in_train])
    
    X_test <- X[in_test, ]
    Y0.hat <- predict(Y0.hat.model, newdata = X_test)$predictions
    Y1.hat <- predict(Y1.hat.model, newdata = X_test)$predictions
    W.hat <- predict(W.hat.model, newdata = X_test)$predictions
    
    W_test <- W[in_test]
    Y.hat <- W_test * Y1.hat + (1 - W_test) * Y0.hat
    
    cate <- Y1.hat - Y0.hat
    po <- cate + ((Y[in_test] - Y.hat) * (W_test - W.hat)) / (W.hat * (1 - W.hat))
    
    list(po = po, cate = cate, Y0.hat = Y0.hat, W.hat = W.hat)
  }) %>% bind_rows()
  
  # tau ----
  tau <- cross_fits$cate
  
  #get out the other things
  W.hat <- cross_fits$W.hat
  Y0.hat <- cross_fits$Y0.hat
  
  # BLP test ----
  BLP_tests <- lapply(seq_len(n_folds), function(fold) {
    in_fold <- fold_indices == fold
    blp_test <- BLP(Y[in_fold], W[in_fold], W.hat[in_fold], Y0.hat[in_fold], tau[in_fold])$coefficients[,c(1,4)]
    return(blp_test)
  })
  # BLP on the whole dataset
  BLP_whole <- BLP(Y, W, W.hat, Y0.hat, tau)$coefficients[,c(1,4)]
  # collate the p-values from the folds
  HTE_pval <- lapply(BLP_tests, function(x) {
    p_val <- x[4,2]
  }) %>% unlist()
  
  # confidence intervals ----
  metaplan <- plan(multicore, workers = workers)
  on.exit(plan(metaplan), add = T)
  draws <- future_map(seq_len(B), function(b) {
    # get your half samples
    half_samples <- lapply(seq_len(n_folds), function(fold) {
      full <- sum(fold_indices == fold)
      half <- c(rep(F, full))
      half[sample(1:full, floor(full/2), replace = F)] <- T
      return(half)
    }) %>% unlist()
    
    # compute the CATEs using the half kept
    tau_half <- lapply(seq_len(n_folds), function(fold) {
      in_train <- half_samples & (fold_indices != fold) # half samples not not in the fold
      train0 <- in_train & W==0 # train control
      train1 <- in_train & W==1 # train treated
      in_fold <- fold_indices == fold
      
      Y0.hat.model <- regression_forest(X[train0, ], Y[train0])
      Y1.hat.model <- regression_forest(X[train1, ], Y[train1])
      
      Y0.hat <- predict(Y0.hat.model, newdata = X[in_fold,])$predictions
      Y1.hat <- predict(Y1.hat.model, newdata = X[in_fold,])$predictions
      
      tau_half_est <- Y1.hat - Y0.hat
      return(tau_half_est)
    }) %>% unlist() %>% unname()
    
    # construct root
    half_root <- tau - tau_half
    return(half_root)
  }, .options = furrr_options(seed = T))
  plan(metaplan)
  
  draws <- do.call(cbind, draws)
  
  # getting the confidence intervals from summary statistics of draws
  lambda_hat <- apply(draws, 1, var)
  normalized <- abs(draws)/(sqrt(lambda_hat))
  col_max    <- apply(normalized, 2, max)
  S_star     <- quantile(col_max, 0.975) # for 95% confidence intervals
  
  # now get the confidence intervals
  res <- data.frame(tau = tau,
                    lb = tau - sqrt(lambda_hat)*S_star,
                    ub = tau + sqrt(lambda_hat)*S_star)
  
  # TE-VIMS ----
  po <- cross_fits$po
  
  te_vims <- NULL
  if (any( HTE_pval < 0.1) | BLP_whole[4,2] < 0.1) {
    covariates <- colnames(X)
    sub_taus <- matrix(nrow = n, ncol = length(covariates))
    colnames(sub_taus) <- covariates
    for (i in seq_along(covariates)) {
      cov <- covariates[i]
      new_X <- as.matrix(X[, -i])
      
      #re estimate cate with new_X
      for (fold in seq_len(n_folds)) {
        in_train <- fold_indices != fold
        in_fold <- !in_train
        
        # fit tau_s by regressing tau on the new_X
        sub_model <- regression_forest(as.matrix(new_X[in_train,]), tau[in_train])
        
        sub_taus[in_fold, i] <- predict(sub_model, newdata = as.matrix(new_X[in_fold,]))$predictions
      }
    }
    
    # Compute TE-VIMs
    ate <- sum(po) / n
    
    r_ate <- (po - ate)^2
    r_tau <- (po - tau)^2
    
    te_vims <- apply(sub_taus, 2, function(sub_tau) {
      r_subtau <- (po - sub_tau)^2
      
      # evaluate TE-VIM (Theta_s in the paper)
      tevim <- sum(r_subtau - r_tau) / n
      
      infl <- r_subtau - r_tau - tevim
      std_err <- sqrt(sum(infl^2)) / n
      
      list(tevim = tevim, std_err = std_err)
    }) %>% simplify2array()
    
    te_vims <- as.data.frame(te_vims)
  }
  
  return(list(tau = res, BLP_tests = BLP_tests, BLP_whole = BLP_whole, te_vims = te_vims, draws = draws))
}
