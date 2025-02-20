######################
# title: model for running the T-learner with a random forest
# date started: 08/01/25
# date finished:
# author: Ellie Van Vogt
#####################

# use the cross fitting code from TE-VIMS package to get the crossfitted CFs...

# Load necessary packages
library(doParallel)
library(dplyr)
library(grf)

# paths ----
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# Reading arguments
args <- commandArgs(trailingOnly = TRUE)
scenario <- as.character(args[1])

# load in the data
datasets <- readRDS(paste0(c("live/data/", scenario,".rds"), collapse = ""))
n_datasets <- length(datasets)

# Set the number of cores for parallel proccess
n_cores <- 10
n_folds <- 10  # Number of folds for cross-fitting

# Run in parallel for all datasets
process_dataset <- function(data) {
  X <- as.matrix(data[, -c(1:4)])  # Keeping only the covariates
  Y <- data$out
  W <- data$trt
  
  n <- dim(X)[1]
  
  # Create fold indices
  fold_indices <- sort(seq(n) %% n_folds) + 1
  fold_list <- unique(fold_indices)
  fold_pairs <- utils::combn(fold_list, 2, simplify = FALSE)  # Fold pairs
  
  # Perform cross-fits
  cross_fits <- lapply(fold_pairs, function(fold_pair) {
    in_train <- !fold_indices %in% fold_pair
    
    Y.hat.model <- regression_forest(X[in_train,], Y[in_train])
    W.hat.model <- regression_forest(X[in_train,], W[in_train])
    
    Y.hat <- predict(Y.hat.model, newdata = X[!in_train,])
    W.hat <- predict(W.hat.model, newdata = X[!in_train,])
    
    list(Y.hat = Y.hat, W.hat = W.hat)
  })
  
  # Collate Y.hat and W.hat matrices
  Y.hat_matrix <- lapply(fold_list, function(fold) {
    Y_hat <- rep(NA, length.out = n)
    for (j in seq_along(fold_pairs)) {
      fold_pair <- fold_pairs[[j]]
      if (fold %in% fold_pair) next
      
      Y_hat[fold_indices %in% fold_pair] <- cross_fits[[j]]$Y.hat[[1]]
    }
    Y_hat[fold_indices == fold] <- NA
    Y_hat
  }) %>% simplify2array()
  
  W.hat_matrix <- lapply(fold_list, function(fold) {
    W_hat <- rep(NA, length.out = n)
    for (j in seq_along(fold_pairs)) {
      fold_pair <- fold_pairs[[j]]
      if (fold %in% fold_pair) next
      
      W_hat[fold_indices %in% fold_pair] <- cross_fits[[j]]$W.hat[[1]]
    }
    W_hat[fold_indices == fold] <- NA
    W_hat
  }) %>% simplify2array()
  
  # Create causal forests
  forest_list <- lapply(fold_list, function(fold) {
    in_train <- fold_indices != fold
    causal_forest(X[in_train,], Y[in_train], W[in_train], Y.hat_matrix[in_train, fold], W.hat_matrix[in_train, fold])
  })
  
  list(W.hat_matrix = W.hat_matrix, Y.hat_matrix = Y.hat_matrix, forest_list = forest_list)
}

t0 <- Sys.time()
# Initialize parallel backend
cl <- makeCluster(n_cores)
registerDoParallel(cl)

results <- foreach(data = datasets, .combine = list, .multicombine = TRUE, .packages = c("grf", "dplyr")) %dopar% {
  process_dataset(data)
}
# Stop parallel backend
stopCluster(cl)

t1 <- Sys.time()
print(t1-t0)
# Save results
saveRDS(results, paste0(c("ephemeral/null/", scenario, "/CF/", "causal_forests.rds"), collapse = ""))
