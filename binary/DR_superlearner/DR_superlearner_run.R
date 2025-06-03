require(SuperLearner, GenericML, dplyr, furrr)
DR_SL_output <- function(data, n_folds, scenario, B, workers, sl_lib) {
  
  X <- data[, -c(1:2)]  # Keeping only the covariates
  Y <- data$Y
  W <- data$W
  
  n <- dim(X)[1]
  
  # Create fold indices
  fold_indices <- sort(seq(n) %% n_folds) + 1
  fold_list <- unique(fold_indices)
  fold_pairs <- utils::combn(fold_list, 2, simplify = FALSE)  # Fold pairs
  
  # Perform cross-fits
  cross_fits <- lapply(fold_pairs, function(fold_pair) {
    in_train <- !(fold_indices %in% fold_pair) #train
    in_test <- which(!in_train)  #test
    X_train <- X[in_train, ]
    X_W_train <- cbind(W = W[in_train], X_train)
    
    Y.hat.model <- SuperLearner(Y = Y[in_train], X = X_W_train, family = binomial(), SL.library = sl_lib, method = method.NNloglik())
    W.hat.model <- SuperLearner(W[in_train], X[in_train, ], family = binomial(), SL.library = sl_lib, method = method.NNloglik())
    X_test <- X[in_test, ]
    X0_test <- cbind(W = 0, X_test)
    X1_test <- cbind(W = 1, X_test)
    Y0.hat <- predict(Y.hat.model, newdata = X0_test)$pred
    Y1.hat <- predict(Y.hat.model, newdata = X1_test)$pred
    W.hat <- predict(W.hat.model, newdata = X_test)$pred
    
    # check if W.hat model failed - then we need to just use the mean as superlearner should be defaulting to the mean!
    if (all(W.hat == 0)) {
      warning("All SuperLearner weights are zero for W.hat. Replacing with constant mean(W).")
      W.hat <- rep(mean(W[in_test]), length = length(in_test))
    }
    
    W_test <- W[in_test]
    Y.hat <- W_test * Y1.hat + (1 - W_test) * Y0.hat
    
    cate <- Y1.hat - Y0.hat
    po <- cate + ((Y[in_test] - Y.hat) * (W_test - W.hat)) / (W.hat * (1 - W.hat))
    
    list(po = po, Y0.hat = Y0.hat, W.hat = W.hat)
  })
  warnings()
  
  # Collate matrices
  po_matrix <- collate_predictions(seq_len(n_folds), fold_pairs, fold_indices, cross_fits, "po")
  Y0.hat_matrix <- collate_predictions(seq_len(n_folds), fold_pairs, fold_indices, cross_fits, "Y0.hat")
  W.hat_matrix <- collate_predictions(seq_len(n_folds), fold_pairs, fold_indices, cross_fits, "W.hat")
  
  # tau ----
  
  tau <- c(rep(NA, n))
  for (fold in seq_along(fold_list)) {
    in_train <- fold_indices != fold
    in_fold <- !in_train
    forest <- SuperLearner(po_matrix[in_train, fold], X[in_train,], family = gaussian(), SL.library = sl_lib)
    tau[in_fold] <- predict(forest, newdata = X[in_fold, ])$pred
  }
  
  # get BLP ----
  Y0.hat <- Y0.hat_matrix %>% rowMeans(,na.rm = T)
  W.hat <- W.hat_matrix %>% rowMeans(, na.rm = T)
  BLP_tests <- lapply(seq_len(n_folds), function(fold) {
    in_fold <- fold_indices == fold
    # BLP is collinear sometimes so the test doesn't work :(
    result <- tryCatch({
      blp_test <- BLP(Y[in_fold], W[in_fold], W.hat[in_fold], Y0.hat[in_fold], tau[in_fold])$coefficients[,c(1,4)]
      return(blp_test)
    },
    error = function(e) {
      # get matrix of the same dim with just NAs to make next steps work
      return(matrix(NA, nrow = 4, ncol = 2))
    })
    
    return(result)
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
    half_samples <- lapply(fold_list, function(fold) {
      full <- sum(fold_indices == fold)
      half <- c(rep(F, full))
      half[sample(1:full, floor(full/2), replace = F)] <- T
      return(half)
    }) %>% unlist()
    
    # compute the CATEs using the half kept
    tau_half <- lapply(fold_list, function(fold) {
      in_train <- half_samples & (fold_indices != fold) # half samples not not in the fold
      in_fold <- fold_indices == fold
      DR_rf <- SuperLearner(po_matrix[in_train, fold], X[in_train,], family = gaussian(), SL.library = sl_lib)
      tau_half_est <- predict(DR_rf, newdata = X[in_fold,])$pred
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
  te_vims <- NULL
  if (any(HTE_pval < 0.1, na.rm = T) | BLP_whole[4,2] < 0.1) {
    print("HTE evidence - running TE-VIMs")
    po <- po_matrix %>% rowMeans(, na.rm = T)
    covariates <- colnames(X)
    sub_taus <- matrix(nrow = n, ncol = length(covariates))
    colnames(sub_taus) <- covariates
    for (i in seq_along(covariates)) {
      cov <- covariates[i]
      new_X <- X %>% select(- cov)
      
      for (fold in seq_along(fold_list)) {
        in_train <- fold_indices != fold
        in_fold <- !in_train
        if (length(covariates) == 2) {
          sl_lib_2 <- sl_lib[!sl_lib %in% "SL.glmnet"]
        } else {
          sl_lib_2 <- sl_lib
        }
        DR_sub <- SuperLearner(po_matrix[in_train, fold], new_X %>% filter(in_train), SL.library = sl_lib_2)
        sub_taus[in_fold, i] <- predict(DR_sub, newdata = new_X %>% filter(in_fold))$pred
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
  } else {
    print("not HTE signal - no TE-VIMs")
  }
  
  return(list(tau = res, BLP_tests = BLP_tests, BLP_whole = BLP_whole, te_vims = te_vims, draws = draws))
}