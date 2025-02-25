######################
# title: TE-VIMs for the DR learner
# date started: 08/01/25
# date finished:
# author: Ellie Van Vogt
#####################
set.seed(1998)
# use the cross fitting code from TE-VIMS package to get the crossfitted CFs...

# Load necessary packages
library(doParallel)
library(dplyr)
library(grf)
library(syrup)

# paths ----
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# functions
source("live/scripts/functions/collate_predictions.R")

# Read args and set params
args <- commandArgs(trailingOnly = TRUE)
scenario <- as.character(args[1])
n <- as.numeric(args[2])
n_cores <- as.numeric(args[3]) 
n_folds <- 10 

# load in the data
datasets <- readRDS(paste0(c("live/data/", scenario, "_", n, ".rds"), collapse = ""))
datasets <- lapply(datasets, `[[`, 1) # just want the data not the truth

# Run in parallel for all datasets
te_vims_DR <- function(data) {
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
    
    list(po = po, Y.hat = Y.hat, W.hat = W.hat)
  })
  
  # Collate cross-fit pseudo outcomes
  po_matrix <- collate_predictions(seq_len(n_folds), fold_pairs, fold_indices, cross_fits, "po")
  
  # out of sample po estimates
  po <- rowMeans(po_matrix, na.rm = T)
  
  # our taus we are comparing against for TE-VIMs
  tau <- c(rep(NA, n))
  for (fold in seq_along(fold_list)) {
    in_train <- fold_indices != fold
    in_fold <- !in_train
    forest <- regression_forest(X[in_train,], po[in_train])
    tau[in_fold] <- predict(forest, newdata = X[in_fold, ])$predictions
  }
  
  # Rerunning the forests with leave-one-out method
  covariates <- colnames(X)
  sub_taus <- matrix(nrow = n, ncol = length(covariates))
  colnames(sub_taus) <- covariates
  for (i in seq_along(covariates)) {
    cov <- covariates[i]
    new_X <- as.matrix(X[, -i])
    
    for (fold in seq_along(fold_list)) {
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
    
  
  return(te_vims)
}


# Parallelise the function
t0 <- Sys.time()
results <- mclapply(datasets, te_vims_DR, mc.cores = n_cores)
t1 <- Sys.time()
print(t1-t0)


# Save results
saveRDS(results, paste0(c("live/results/", scenario, "/", n,  "/DR_RF/", "DR_rf_tevims.RDS"), collapse = ""))
