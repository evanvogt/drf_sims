###################
# Required packages
###################
require(coin)
require(grf)
require(future)
require(furrr)
require(dplyr)

###################
# Main function for running causal forests - sub-distribution and cause specific approaches
###################

run_all_cate_methods_survival <- function(data, n_folds = 10, workers = 5, 
                                          horizon = NULL) {
  # Set up parallel processing plan
  metaplan <- plan(multisession, workers = workers)
  on.exit(plan(metaplan), add = TRUE)
  # Data preparation
  X <- as.matrix(data[, !(names(data) %in% c("Y", "D", "W"))])
  Y <- data$Y  # Event time
  D <- data$D  # Event type (0 = censored, 1 = event of interest, 2 = competing event)
  W <- data$W  # Treatment indicator
  n_obs <- nrow(X)
  # If no horizon specified, use maximum event time
  if (is.null(horizon)) {
    horizon <- max(Y[D != 0])  # Use max observed event time
    warning("No horizon specified. Using maximum observed event time: ", horizon)
  }
  # Create fold indices 
  fold_indices <- sort(seq(n_obs) %% n_folds) + 1
  fold_list <- unique(fold_indices)
  fold_pairs <- utils::combn(fold_list, 2, simplify = FALSE)
  # Results container
  results <- list()
  # 1. CSF - subdistribution weighting approach ----
  cat("Running CSF with subdistribution weights...\n")
  results$causal_survival_forest_sub_dist <- run_csf_sub_dist(
    X, Y, D, W, fold_indices, fold_list, n_folds, horizon
  )
  # 2. SF - subdistribution weighting approach ----
  cat("Running CSF with cause-specific weighting approach...\n")
  results$causal_survival_forest_cause_spec <- run_csf_cause_spec(
    X, Y, D, W, fold_indices, fold_list, n_folds, horizon
  )
  return(results)
}

###################
# Helper functions - censoring
###################
handle_censoring_sub_dist <- function(X, Y, D, W, horizon, in_train) {
  n_obs <- nrow(X)
  # Create pseudo-outcomes for competing risks (Aalen-Johansen style)
  Y_pseudo <- Y
  Y_pseudo[D == 2] <- horizon + 1  # Move competing events beyond horizon
  D_pseudo <- as.integer(D == 1)   # Binary indicator for EOI or censor
  # Check for censoring
  has_censoring <- any(D[in_train] == 0 & Y[in_train] < horizon)
  if (!has_censoring) {
    # No censoring during follow-up - simple case
    return(list(
      Y_pseudo = Y_pseudo,
      D_pseudo = D_pseudo,
      sample_weights = rep(1, n_obs),
      observed_events = rep(TRUE, n_obs),
      has_censoring = FALSE
    ))
  } else {
    # Handle censoring using IPCW
    X_train <- X[in_train,]
    W_train <- W[in_train]
    Y_train <- Y[in_train]
    D_train <- D[in_train]
    
    # Fit time-to-censoring distribution
    sf_censor <- survival_forest(cbind(X_train, W_train), 
                                 Y = Y_train, 
                                 D = (D_train == 0))
    # Predict censoring probabilities for all observations (removal of in-fold predictions happens outside of function)
    sf_pred <- predict(sf_censor, 
                       newdata = cbind(X, W),
                       failure.times = pmin(Y, horizon), 
                       prediction.times = "time")
    censoring_prob <- sf_pred$predictions
    # Identify observed events (not censored before horizon)
    observed_events <- (D != 0 | (D == 0 & Y >= horizon))
    # Calculate IPCW weights
    sample_weights <- rep(1, n_obs)
    sample_weights[observed_events] <- 1 / pmax(censoring_prob[observed_events],1e-3)
    return(list(
      Y_pseudo = Y_pseudo,
      D_pseudo = D_pseudo,
      sample_weights = sample_weights,
      observed_events = observed_events,
      has_censoring = TRUE,
      censoring_prob = censoring_prob
    ))
  }
}

handle_censoring_cause_spec <- function(X, Y, D, W, horizon, in_train) {
  n_obs <- length(Y)
  D_event <- as.integer(D == 1 & Y <= horizon)
  has_censoring <- any(D[in_train] == 0 & Y[in_train] <= horizon)
  has_competing <- any(D[in_train] == 2 & Y[in_train] <= horizon)
  X_train <- X[in_train,]
  W_train <- W[in_train]
  Y_train <- Y[in_train]
  D_train <- D[in_train]
  sample_weights <- rep(1, n_obs)
  censor_surv_prob <- rep(1, n_obs)
  comp_surv_prob <- rep(1, n_obs)
  
  if (has_censoring) {
    sf_censor <- survival_forest(cbind(X_train, W_train),
                                 Y = Y_train,
                                 D = as.integer(D_train == 0))
    censor_surv_prob <- predict(sf_censor,
                                newdata = cbind(X, W), 
                                failure.times = pmin(Y, horizon),
                                prediction.times = "time")$predictions
    censor_surv_prob <- pmax(censor_surv_prob, 1e-3)
  } 
  if (has_competing) {
    sf_comp <- survival_forest(cbind(X_train, W_train),
      Y = Y_train,
      D = as.integer(D_train == 2)
    )
    comp_surv_prob <- predict(sf_comp, 
                              newdata = cbind(X, W),
                              prediction.times = "time",
                              failure.times = pmin(Y, horizon))$predictions
    comp_surv_prob <- pmax(comp_surv_prob, 1e-3)
  }
  combined_surv <- pmax(censor_surv_prob * comp_surv_prob, 1e-3)
  observed_at_time <- (D != 0) | (Y >= horizon)
  sample_weights[observed_at_time] <- 1 / combined_surv[observed_at_time]
  return(list(
    Y_event = Y,
    D_event = D_event,
    sample_weights = sample_weights,
    observed_at_horizon = observed_at_time,
    has_censoring = has_censoring,
    has_competing = has_competing
  ))
}

###################
# Method 1: CSF A-J / sub-dist
###################
run_csf_sub_dist <- function(X, Y, D, W, fold_indices, fold_list, n_folds, horizon) {
  n_obs <- nrow(X)
  # Cross-validated tau estimation with fold-specific censoring weights
  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold <- !in_train
    # Handle censoring 
    censoring_info <- handle_censoring_sub_dist(X, Y, D, W, horizon, in_train)
    Y_pseudo <- censoring_info$Y_pseudo
    D_pseudo <- censoring_info$D_pseudo
    sample_weights <- censoring_info$sample_weights
    observed_events <- censoring_info$observed_events
    if (censoring_info$has_censoring) {
      # With censoring - use only observed events for training
      train_obs <- in_train & observed_events
      if (sum(train_obs) < 20) {
        return(list(fold = fold, tau = rep(NA, sum(in_fold))))
      }
      forest <- causal_survival_forest(
        X = X[train_obs,], 
        Y = Y_pseudo[train_obs], 
        D = D_pseudo[train_obs], 
        W = W[train_obs],
        sample.weights = sample_weights[train_obs],
        horizon = horizon
      )
    } else {
      # No censoring - use all training observations
      forest <- causal_survival_forest(
        X = X[in_train,], 
        Y = Y_pseudo[in_train], 
        D = D_pseudo[in_train], 
        W = W[in_train],
        horizon = horizon
      )
    }
    # Predict for test fold
    tau_fold <- predict(forest, newdata = X[in_fold,])$predictions
    list(fold = fold, tau = tau_fold)
  }, .options = furrr_options(seed = TRUE))
  # Reconstruct tau vector
  tau <- rep(NA, n_obs)
  for (result in tau_results) {
    fold <- result$fold
    in_fold <- fold_indices == fold
    tau[in_fold] <- result$tau
  }
  return(list(
    tau = data.frame(tau = tau),
    method = "causal_survival_forest_sub_dist"
  ))
}

###################
# Method 2: CSF F-G / cause-spec
###################
run_csf_cause_spec <- function(X, Y, D, W, fold_indices, fold_list, n_folds, horizon) {
  n_obs <- nrow(X)
  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold <- !in_train
    cr_info <- handle_censoring_cause_spec(X, Y, D, W, horizon, in_train)
    Y_event <- cr_info$Y_event
    D_event <- cr_info$D_event
    sample_weights <- cr_info$sample_weights
    if (cr_info$has_censoring || cr_info$has_competing) {
      train_obs <- in_train & cr_info$observed_at_horizon
      if (sum(train_obs) < 20) {
        return(list(fold=fold, tau=rep(NA,sum(in_fold))))
      }
      forest <- causal_survival_forest(X = X[train_obs,],
                                       Y = Y_event[train_obs],
                                       D = D_event[train_obs],
                                       W = W[train_obs],
                                       sample.weights = sample_weights[train_obs],
                                       horizon = horizon)
      } else {
      forest <- causal_survival_forest(X = X[in_train,],
                                       Y=Y_event[in_train],
                                       D=D_event[in_train],
                                       W=W[in_train],
                                       horizon=horizon)
    }
    tau_fold <- predict(forest, newdata=X[in_fold,,drop=FALSE])$predictions
    list(fold=fold, tau=tau_fold)
  }, .options=furrr_options(seed=TRUE))
  tau <- rep(NA,n_obs)
  for (res in tau_results) {
    tau[fold_indices==res$fold] <- res$tau
  }
  return(list(tau=data.frame(tau=tau),
              method="causal_survival_forest_cause_spec"))
}
