######################
# title: TE-vims for the T-learner
# date started: 08/01/25
# date finished:
# author: Ellie Van Vogt
#####################
set.seed(1998)
# adapting O hines github code

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

# Read args and set parameters
args <- commandArgs(trailingOnly = TRUE)
scenario <- as.character(args[1])
n <- as.numeric(args[2])
n_cores <- as.numeric(args[3])
n_folds <- 10 

# load in the data
datasets <- readRDS(paste0(c("live/data/", scenario, "_", n, ".rds"), collapse = ""))
datasets <- lapply(datasets, `[[`, 1) # just want the data not the truth



# TE-VIMs for a single dataset
T_RF_tevim <- function(data) {
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
    
    list(po = po, cate = cate)
  })
  
  # Collate po and tau
  po <- rep(NA, length.out = n)
  tau <- rep(NA, length.out = n)
  for (fold in seq_len(n_folds)) {
    in_fold <- fold_indices == fold
    po[in_fold] <- cross_fits[[fold]]$po
    tau[in_fold] <- cross_fits[[fold]]$cate
  }
  
  # getting the subtaus for removing each variable
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
      train0 <- in_train & W==0 # train control
      train1 <- in_train & W==1 # train treated
      
      # re train Y models on new X
      Y0.hat.model <- regression_forest(as.matrix(new_X[train0, ]), Y[train0])
      Y1.hat.model <- regression_forest(as.matrix(new_X[train1, ]), Y[train1])

      Y0.hat <- predict(Y0.hat.model, newdata = as.matrix(new_X[in_fold, ]))$predictions
      Y1.hat <- predict(Y1.hat.model, newdata = as.matrix(new_X[in_fold, ]))$predictions
      
      # cate in fold
      cate <- Y1.hat - Y0.hat
      
      sub_taus[in_fold, i] <- cate
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
results <- mclapply(datasets, T_RF_tevim, mc.cores = n_cores)
t1 <- Sys.time()
print(t1-t0)

# Save results
saveRDS(results, paste0(c("live/results/", scenario, "/", n, "/T_RF/", "T_RF_te_vims.RDS"), collapse = ""))