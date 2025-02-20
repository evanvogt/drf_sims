######################
# title: model for running the DR learner with cross fitting
# date started: 08/01/25
# date finished:
# author: Ellie Van Vogt
#####################

# running the models mainly to check running times

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

# Set the number of cores for parallel proccess
n_cores <- 40
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
    
    Y.hat.model <- regression_forest(cbind(W[in_train], X[in_train,]), Y[in_train])
    W.hat.model <- regression_forest(X[in_train,], W[in_train])
    
    Y0.hat <- predict(Y.hat.model, newdata = cbind(W = 0, X[!in_train,]))
    Y1.hat <- predict(Y.hat.model, newdata = cbind(W = 1, X[!in_train,]))
    W.hat <- predict(W.hat.model, newdata = X[!in_train,])
    
    Y.hat <- if_else(W[!in_train] > 0.5, Y1.hat, Y0.hat)
    cate <- Y1.hat - Y0.hat
    po <- cate + (Y[!in_train] - Y.hat) * (W[!in_train] - W.hat) / (W.hat *(1 - W.hat))
    
    list(po = po, Y.hat = Y.hat, W.hat = W.hat)
  })
  
  # Collate cross-fit pseudo outcomes
  po_matrix <- lapply(fold_list, function(fold) {
    po_hat <- rep(NA, length.out = n)
    for (i in seq_along(fold_pairs)) {
      fold_pair <- fold_pairs[[i]]
      if (fold %in% fold_pair) {
        next
      }
      po_hat[fold_indices %in% fold_pair] <- cross_fits[[i]]$po[[1]]
    }
    po_hat[fold_indices == fold] <- NA
    po_hat
  }) %>% simplify2array()
  
  # out of sample po estimates
  po <- rowMeans(po_matrix, na.rm = T)
  
  # last stage regression to get the taus

  DR_list <- lapply(fold_list, function(fold) {
    in_train <- fold_indices != fold
    in_fold <- !in_train
    DR_rf <- regression_forest(X[in_train,], po[in_train])
    DR_rf
  })
  
  list(po_matrix = po_matrix, DR_list = DR_list)
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
saveRDS(results, paste0(c("ephemeral/null/", scenario, "/DR learner/", "DR_model.rds"), collapse = ""))
