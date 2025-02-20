######################
# title: model for running the causal forest with cross-fitting
# date started: 08/01/25
# date finished:
# author: Ellie Van Vogt
#####################
set.seed(1998)
# getting out the CATEs and the BLPs

# Load necessary packages
library(foreach)
library(doParallel)
library(purrr)
library(grf)

# paths ----
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

# Set the number of cores for parallel proccess
n_cores <- 10 #floor(future::availableCores() *0.9)
n_folds <- 10  # Number of folds for cross-fitting

# Run in parallel for all datasets
cate_blp <- function(data) {
  X <- as.matrix(data[, -c(1:2)])  # Keeping only the covariates
  Y <- data$Y
  W <- data$W
  
  n <- dim(X)[1]
  
  # Create fold indices
  fold_indices <- sort(seq(n) %% n_folds) + 1
  fold_pairs <- combn(seq_len(n_folds), 2, simplify = FALSE)
  
  # Perform cross-fits
  cross_fits <- lapply(fold_pairs, function(fold_pair) {
    in_train <- !fold_indices %in% fold_pair
    
    Y.hat.model <- regression_forest(X[in_train,], Y[in_train])
    W.hat.model <- regression_forest(X[in_train,], W[in_train])
    
    Y.hat <- predict(Y.hat.model, newdata = X[!in_train,])$predictions
    W.hat <- predict(W.hat.model, newdata = X[!in_train,])$predictions
    
    list(Y.hat = Y.hat, W.hat = W.hat)
  })

  # Collate Y.hat and W.hat matrices
  Y.hat_matrix <- collate_predictions(seq_len(n_folds), fold_pairs, fold_indices, cross_fits, "Y.hat")
  W.hat_matrix <- collate_predictions(seq_len(n_folds), fold_pairs, fold_indices, cross_fits, "W.hat")
  
  # Create causal forests
  tau <- matrix(NA, n, 3)
  BLP <- numeric(n_folds)
  for (fold in seq_len(n_folds)) {
    in_train <- fold_indices != fold
    in_fold <- !in_train
    
    forest <- causal_forest(X[in_train,], Y[in_train], W[in_train], 
                            Y.hat_matrix[in_train, fold], W.hat_matrix[in_train, fold])
    
    pred <- predict(forest, newdata = X[in_fold,], estimate.variance = TRUE)
    
    tau_est <- cbind(
      pred$predictions,
      pred$predictions - qnorm(0.975) * sqrt(pred$variance.estimates),
      pred$predictions + qnorm(0.975) * sqrt(pred$variance.estimates)
    )
    tau[in_fold, ] <- tau_est
    BLP[fold] <- test_calibration(forest)[2, 4]
  }
  list(tau = tau, BLP = BLP)
}


# Parallelise the function
t0 <- Sys.time()
results <- mclapply(datasets, cate_blp, mc.cores = n_cores)
t1 <- Sys.time()
print(t1-t0)


# Separate results into CATEs and BLP tests
CATEs <- sapply(results, `[[`, "tau")  # Extract the "CATE" element from each result
BLP_tests <- sapply(results, `[[`, "BLP")  # Extract the "BLP_test" element from each result

CATEs <- t(CATEs)
BLP_tests <- t(BLP_tests)

# save output
saveRDS(CATEs, paste0(c("live/results/", scenario, "/", n, "/CF/", "causal_forest_cates.RDS"), collapse = ""))
saveRDS(BLP_tests, paste0(c("live/results/", scenario, "/", n, "/CF/", "causal_forest_BLP.RDS"), collapse = ""))