######################
# title: confidence intervals for the T-learner
# date started: 08/01/25
# date finished:
# author: Ellie Van Vogt
#####################
set.seed(1998)
# simultaneous inference code - my version

# libraries
library(furrr)
library(dplyr)
library(grf)
library(GenericML)

# paths 
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# functions
source("live/scripts/functions/collate_predictions.R")

# Read args and set parameters
args <- commandArgs(trailingOnly = TRUE)
scenario <- as.character(args[1])
n <- as.numeric(args[2])
n_cores <- as.numeric(args[3])
n_folds <- 10 
B <- 200

oldplan <- plan(multisession, workers = n_cores)

# load in the data
datasets <- readRDS(paste0(c("live/data/", scenario, "_", n, ".RDS"), collapse = ""))
datasets <- lapply(datasets, `[[`, 1) # just want the data not the truth


# define a function to get CIs and taus for a single dataset
T_RF_cis <- function(data) {
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
  
  # get the point taus
  tau <- rep(NA, length.out = n)
  for (fold in seq_len(n_folds)) {
    in_fold <- fold_indices == fold
    tau[in_fold] <- cross_fits$cate[in_fold]
  }
  
  # get the BLP
  BLP_tests <- lapply(seq_len(n_folds), function(fold) {
    in_fold <- fold_indices == fold
    blp_test <- BLP(Y[in_fold], W[in_fold], cross_fits$W.hat[in_fold], cross_fits$Y0.hat[in_fold], cross_fits$cate[in_fold])$coefficients[,c(1,4)]
    return(blp_test)
  })
  

  
  t2 <- Sys.time()
  # B half bootstraps - parallelise this to hopefully speed it up
  metaplan <- plan(multicore, workers = 5)
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
  
  t3 <- Sys.time()
  print(t3-t2)
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
  return(list(res= res, draws = draws, BLP = BLP_tests))
}



# Parallelise the function
t0 <- Sys.time()
results <- future_map(datasets, T_RF_cis, .options = furrr_options(seed = T))
t1 <- Sys.time()
print(t1-t0)
plan(oldplan)

#separate results
CATEs <- lapply(results, `[[`, "res") 
BLP_tests <- lapply(results, `[[`, "BLP_tests") 
draws <- lapply(results, `[[`, "draws") 

# Save results
saveRDS(CATEs, paste0(c("live/results/", scenario, "/", n, "/T_RF/", "tau_cis.RDS"), collapse = ""))
saveRDS(BLP_tests, paste0(c("live/results/", scenario, "/", n, "/T_RF/", "BLP.RDS"), collapse = ""))
saveRDS(draws, paste0(c("live/results/", scenario, "/", n, "/T_RF/", "draws.RDS"), collapse = ""))