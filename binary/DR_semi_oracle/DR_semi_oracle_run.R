require(grf, GenericML, dplyr, furrr)
DR_semi_oracle_output <- function(data, n_folds, scenario, B, workers) {

  X <- as.matrix(data[, -c(1:2)])  # Keeping only the covariates
  Y <- data$Y
  W <- data$W
  n <- dim(X)[1]
  
  # specify propensity as 0.5
  W.hat <- rep(0.5, n)
  
  # need to learn Y.hat only - simpler crossfitting
  
  # Create fold indices
  fold_indices <- sort(seq(n) %% n_folds) + 1
  
  # Perform cross-fits
  cross_fits <- lapply(seq_len(n_folds), function(fold) {
    in_train <- !(fold_indices == fold) #train
    in_test <- which(!in_train)  #test
    
    Y.hat.model <- regression_forest(cbind(W[in_train], X[in_train, ]), Y[in_train])
    
    X_test <- X[in_test, ]
    Y0.hat <- predict(Y.hat.model, newdata = cbind(W = 0, X_test))$predictions
    Y1.hat <- predict(Y.hat.model, newdata = cbind(W = 1, X_test))$predictions
    
    W_test <- W[in_test]
    Y.hat <- W_test * Y1.hat + (1 - W_test) * Y0.hat
    W.hat <- rep(0.5, length(in_test))
    cate <- Y1.hat - Y0.hat
    po <- cate + ((Y[in_test] - Y.hat) * (W_test - W.hat)) / (W.hat * (1 - W.hat))
    
    list(po = po, Y0.hat = Y0.hat)
  }) %>% bind_rows()
  
  
  #po and Y0
  po <- cross_fits$po
  Y0.hat <- cross_fits$Y0.hat

  # taus ----
  tau <- c(rep(NA, n))
  for (fold in seq_len(n_folds)) {
    in_train <- fold_indices != fold
    in_fold <- !in_train
    forest <- regression_forest(X[in_train,], po[in_train])
    tau[in_fold] <- predict(forest, newdata = X[in_fold, ])$predictions
  }
  
  # BLP tests ----
  BLP_tests <- lapply(seq_len(n_folds), function(fold) {
    in_fold <- fold_indices == fold
    result <- tryCatch({
      blp_test <- BLP(Y[in_fold], W[in_fold], W.hat[in_fold], Y0.hat[in_fold], tau[in_fold])$coefficients[,c(1,4)]
      return(blp_test)
    },
    error = function(e) {
      # get matrix of the same dim with just NAs to make next steps work
      return(matrix(1, nrow = 4, ncol = 2))
    })
    
    return(result)
  })
  
  # BLP on the whole dataset
  BLP_whole <- BLP(Y, W, W.hat, Y0.hat, tau)$coefficients[,c(1,4)]
  # collate the p-values from the folds
  HTE_pval <- lapply(BLP_tests, function(x) {
    p_val <- x[4,2]
    if(is.na(p_val)) {
      p_val <- 1
    }
    return(p_val)
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
      in_fold <- fold_indices == fold
      DR_rf <- regression_forest(X[in_train,], po[in_train])
      tau_half_est <- predict(DR_rf, newdata = X[in_fold,])
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
  
  
  
  # TE-VIMS -----
  te_vims <- NULL
  if (any( HTE_pval < 0.1) | BLP_whole[4,2] < 0.1) {
    covariates <- colnames(X)
    sub_taus <- matrix(nrow = n, ncol = length(covariates))
    colnames(sub_taus) <- covariates
    for (i in seq_along(covariates)) {
      cov <- covariates[i]
      new_X <- as.matrix(X[, -i])
      
      for (fold in seq_len(n_folds)) {
        in_train <- fold_indices != fold
        in_fold <- !in_train
        DR_sub <- regression_forest(as.matrix(new_X[in_train, ]), po[in_train])
        sub_taus[in_fold, i] <- predict(DR_sub, newdata = as.matrix(new_X[in_fold, ]))$predictions
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