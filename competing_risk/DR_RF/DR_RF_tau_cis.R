######################
# title: getting tau estimates and confidence intervals for the DR learner with RFs
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
source("live/scripts/competing_risk/functions/collate_predictions.R")

# Reading arguments and set params
args <- commandArgs(trailingOnly = TRUE)
scenario <- as.character(args[1])
n <- as.numeric(args[2])
n_cores <- as.numeric(args[3]) 
n_folds <- 10
B <- 200 # number of bootstraps

oldplan <- plan(multisession, workers = n_cores)

# load in the data
datasets <- readRDS(paste0(c("live/data/competing_risk//", scenario, "_", n, ".RDS"), collapse = ""))
truths <- lapply(datasets, `[[`, 2)
datasets <- lapply(datasets, `[[`, 1) # just want the data not the truth

# define a function to get CIs and taus for a single dataset
taus_and_cis <- function(data) {
  X <- as.matrix(data[, -c(1:2)])  # Keeping only the covariates
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
    
    list(po = po, Y0.hat = Y0.hat, W.hat = W.hat)
  })
  
  # Collate matrices
  po_matrix <- collate_predictions(seq_len(n_folds), fold_pairs, fold_indices, cross_fits, "po")
  Y0.hat_matrix <- collate_predictions(seq_len(n_folds), fold_pairs, fold_indices, cross_fits, "Y0.hat")
  W.hat_matrix <- collate_predictions(seq_len(n_folds), fold_pairs, fold_indices, cross_fits, "W.hat")
  
  po <- po_matrix %>% rowMeans(, na.rm = T)
  Y0.hat <- Y0.hat_matrix %>% rowMeans(, na.rm = T)
  W.hat <- W.hat_matrix %>% rowMeans(, na.rm = T)

  # point estimate for tau
  tau <- c(rep(NA, n))
  for (fold in seq_along(fold_list)) {
    in_train <- fold_indices != fold
    in_fold <- !in_train
    forest <- regression_forest(X[in_train,], po[in_train])
    tau[in_fold] <- predict(forest, newdata = X[in_fold, ])$predictions
  }
  
  # get BLP
  BLP_tests <- lapply(seq_len(n_folds), function(fold) {
    in_fold <- fold_indices == fold
    blp_test <- BLP(Y[in_fold], W[in_fold], W.hat[in_fold], Y0.hat[in_fold], tau[in_fold])$coefficients[,c(1,4)]
    return(blp_test)
  })
  
  t2 <- Sys.time()
  # B half bootstraps - parallelise this to hopefully speed it up
  metaplan <- plan(multicore, workers = 5)
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
      DR_rf <- regression_forest(X[in_train,], po[in_train])
      tau_half_est <- predict(DR_rf, newdata = X[in_fold,])
      return(tau_half_est)
    }) %>% unlist() %>% unname()
    
    # construct root
    half_root <- tau - tau_half
    return(half_root)
  }, .options = furrr_options(seed = T))
  plan(metaplan)
  t3 <- Sys.time()
  print(t3-t2)
  
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
results <- future_map(datasets, taus_and_cis, .options = furrr_options(seed = T))
t1 <- Sys.time()
print(t1-t0)
plan(oldplan)


# Separate results
CATEs <- lapply(results, `[[`, "res") 
BLP_tests <- lapply(results, `[[`, "BLP_tests") 
draws <- lapply(results, `[[`, "draws") 

# Save results
saveRDS(CATEs, paste0(c("live/results/competing_risk/", scenario, "/", n, "/DR_RF/", "taus_cis.RDS"), collapse = ""))
saveRDS(BLP_tests, paste0(c("live/results/competing_risk/", scenario, "/", n, "/DR_RF/", "BLP.RDS"), collapse = ""))
saveRDS(draws, paste0(c("live/results/competing_risk/", scenario, "/", n, "/DR_RF/", "draws.RDS"), collapse = ""))