######################
# title: model-free variable importance metric for the causal forest
# date started: 08/01/25
# date finished:
# author: Ellie Van Vogt
#####################
set.seed(1998)
# Libraries
library(doParallel)
library(dplyr)
library(grf)
library(syrup)

# Paths
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


# Parameters to set
n_cores <- 10  
n_folds <- 10  


te_vims_CF <- function(data) {
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
  
  
  # Collate Y.hat and W.hat matrices
  collate_predictions <- function(fold_list, fold_pairs, fold_indices, cross_fits, target) {
    lapply(fold_list, function(fold) {
      predictions <- rep(NA, length.out = length(fold_indices))
      for (j in seq_along(fold_pairs)) {
        fold_pair <- fold_pairs[[j]]
        if (fold %in% fold_pair) next
        
        predictions[fold_indices %in% fold_pair] <- cross_fits[[j]][[target]]
      }
      predictions[fold_indices == fold] <- NA
      predictions
    }) %>% simplify2array()
  }
  
  po_matrix <- collate_predictions(seq_len(n_folds), fold_pairs, fold_indices, cross_fits, "po")
  Y.hat_matrix <- collate_predictions(seq_len(n_folds), fold_pairs, fold_indices, cross_fits, "Y.hat")
  W.hat_matrix <- collate_predictions(seq_len(n_folds), fold_pairs, fold_indices, cross_fits, "W.hat")
  

  tau <- c(rep(NA, n))
  for (fold in seq_along(fold_list)) {
    in_train <- fold_indices != fold
    in_fold <- !in_train
    forest <- causal_forest(X[in_train, ], Y[in_train], W[in_train], Y.hat_matrix[in_train, fold], W.hat_matrix[in_train, fold])
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
      new_forest <- causal_forest(as.matrix(new_X[in_train, ]), Y[in_train], W[in_train], Y.hat_matrix[in_train, fold], W.hat_matrix[in_train, fold])
      sub_taus[in_fold, i] <- predict(new_forest, newdata = as.matrix(new_X[in_fold, ]))$predictions
    }
  }
  
  # Compute TE-VIMs
  po <- rowMeans(po_matrix, na.rm = T)
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
  })
  
  return(te_vims)
}

# Parallelise the function
t0 <- Sys.time()
results <- mclapply(datasets, te_vims_CF, mc.cores = n_cores)
t1 <- Sys.time()
print(t1-t0)

# Save results
saveRDS(results, file = paste0(c("live/results/", scenario, "/", n, "/CF/", "te_vims_list.rds"), collapse = ""))

