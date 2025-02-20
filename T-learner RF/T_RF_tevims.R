######################
# title: TE-vims for the T-learner
# date started: 08/01/25
# date finished:
# author: Ellie Van Vogt
#####################
set.seed(1998)
# simultaneous inference code - my version

# libraries
library(foreach)
library(doParallel)
library(dplyr)
library(grf)
library(syrup)

# paths 
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# functions
source("live/scripts/functions/collate_predictions.R")

# Reading arguments
args <- commandArgs(trailingOnly = TRUE)
scenario <- as.character(args[1])
n <- as.numeric(args[2])

# load in the data
datasets <- readRDS(paste0(c("live/data/", scenario, "_", n, ".rds"), collapse = ""))
datasets <- lapply(datasets, `[[`, 1) # just want the data not the truth

# parameters
n_cores <- 10 #floor(future::availableCores() *0.9)
n_folds <- 10 
B <- 500 # number of bootstraps

# define a function to get CIs and taus for a single dataset
T_RF_cis <- function(data) {
  X <- as.matrix(data[, -c(1:2)])  # Keeping only the covariates
  Y <- data$Y
  W <- data$W
  
  n <- dim(X)[1]
  
  # Create fold indices
  fold_indices <- sort(seq(n) %% n_folds) + 1
  fold_list <- unique(fold_indices)
  fold_pairs <- utils::combn(fold_list, 2, simplify = FALSE)  # Fold pairs
  
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
    
    list(po = po, cate = cate, Y.hat = Y.hat, W.hat = W.hat)
  })
  
  # Collate po and tau
  po <- rep(NA, length.out = n)
  tau <- rep(NA, length.out = n)
  for (fold in seq_len(n_folds)) {
    in_fold <- fold_indices == fold
    po[in_fold] <- cross_fits[[fold]]$po
    tau[in_fold] <- cross_fits[[fold]]$cate
  }
  
  
  
  t2 <- Sys.time()
  # B half bootstraps
  draws <- replicate(B, {
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
      DR_rf <- regression_forest(X[in_train,], po[in_train])
      tau_half_est <- predict(DR_rf, newdata = X[in_fold,])
      return(tau_half_est)
    }) %>% unlist() %>% unname()
    
    # construct root
    half_root <- tau - tau_half
    return(half_root)
  })
  t3 <- Sys.time()
  print(t3-t2)
  
  # getting the confidence intervals from summary statistics of draws
  lambda_hat <- apply(draws, 1, var)
  normalized <- abs(draws)/(sqrt(lambda_hat))
  col_max    <- apply(normalized, 2, max)
  S_star     <- quantile(col_max, 0.975) # for 95% confidence intervals
  
  # now get the confidence intervals
  res <- data.frame(tau = tau,
                    lb = tau - sqrt(lambda_hat)*S_star,
                    ub = tau + sqrt(lambda_hat)*S_star)
  return(res)
}



# Parallelise the function
t0 <- Sys.time()
results <- mclapply(datasets, taus_and_cis, mc.cores = n_cores)
t1 <- Sys.time()
print(t1-t0)

# Save results
saveRDS(results, paste0(c("live/results/", scenario, "/", n, "/DR learner/", "DR_rf_taus_cis.RDS"), collapse = ""))