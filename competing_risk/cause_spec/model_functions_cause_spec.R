###################
# Required packages
###################
require(grf)
require(coin)
require(future)
require(furrr)
require(dplyr)

###################
# Main function for running CATE models on competing risks survival data
###################
run_all_cate_methods_survival <- function(data, n_folds = 10, workers = 5, 
                                          horizon = NULL, fmla_info = NULL) {
  metaplan <- plan(multisession, workers = workers)
  on.exit(plan(metaplan), add = TRUE)
  X <- as.matrix(data[, !(names(data) %in% c("Y", "D", "W"))])
  Y <- data$Y
  D <- data$D
  W <- data$W
  n_obs <- nrow(X)
  if (is.null(horizon)) {
    horizon <- max(Y[D != 0])
    warning("No horizon specified. Using maximum observed event time: ", horizon)
  }
  fold_indices <- sort(seq(n_obs) %% n_folds) + 1
  fold_list <- unique(fold_indices)
  fold_pairs <- utils::combn(fold_list, 2, simplify = FALSE)
  results <- list()
  cat("Running Causal Survival Forest...\n")
  results$causal_survival_forest <- run_causal_survival_forest(
    X, Y, D, W, fold_indices, fold_list, n_folds, horizon
  )
  if (!is.null(fmla_info)) {
    cat("Running DR Oracle...\n")
    results$dr_oracle <- run_dr_oracle_survival(
      X, Y, D, W, fmla_info, fold_indices, fold_list, n_folds, horizon
    )
  }
  cat("Running DR Semi-Oracle...\n")
  results$dr_semi_oracle <- run_dr_semi_oracle_survival(
    X, Y, D, W, fold_indices, fold_list, n_folds, horizon
  )
  cat("Running DR Survival Forest (grf only)...\n")
  results$dr_survival_rf <- run_dr_survival_grf(
    X, Y, D, W, fold_indices, fold_list, n_folds, horizon
  )
  return(results)
}

###################
# Helper for IPCW weighting for censoring and competing risk
###################
handle_censoring_and_competing_risks <- function(X, Y, D, W, horizon, in_train) {
  n_obs <- length(Y)
  D_event <- as.integer(D == 1 & Y <= horizon)
  has_censoring <- any(D[in_train] == 0 & Y[in_train] <= horizon)
  has_competing <- any(D[in_train] == 2 & Y[in_train] <= horizon)
  X_train <- X[in_train, , drop=FALSE]
  W_train <- W[in_train]
  Y_train <- Y[in_train]
  D_train <- D[in_train]
  sample_weights <- rep(1, n_obs)
  censor_surv_prob <- rep(1, n_obs)
  comp_surv_prob <- rep(1, n_obs)
  
  if (has_censoring) {
    sf_censor <- survival_forest(
      cbind(X_train, W_train),
      Y = Y_train,
      D = as.integer(D_train == 0)
    )
    eval_times <- pmin(Y, horizon)
    censor_surv_prob <- predict(sf_censor, 
                                newdata = cbind(X, W), prediction.times = "time",
                                failure.times = eval_times)$predictions
    censor_surv_prob <- pmax(censor_surv_prob, 1e-6)
  } 
  if (has_competing) {
    sf_comp <- survival_forest(
      cbind(X_train, W_train),
      Y = Y_train,
      D = as.integer(D_train == 2)
    )
    eval_times <- pmin(Y, horizon)
    comp_surv_prob <- predict(sf_comp, 
                              newdata = cbind(X, W),predict.times = "time",
                              failure.times = eval_times)$predictions
    comp_surv_prob <- pmax(comp_surv_prob, 1e-6)
  }
  combined_surv <- censor_surv_prob * comp_surv_prob
  observed_at_time <- (D != 0) | (Y >= horizon)
  sample_weights[observed_at_time] <- 1 / combined_surv[observed_at_time]
  sample_weights <- pmin(sample_weights, 10)
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
# Method 1: Causal Survival Forest 
###################
run_causal_survival_forest <- function(X, Y, D, W, fold_indices, fold_list, n_folds, horizon) {
  n_obs <- nrow(X)
  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold <- !in_train
    cr_info <- handle_censoring_and_competing_risks(X, Y, D, W, horizon, in_train)
    Y_event <- cr_info$Y_event
    D_event <- cr_info$D_event
    sample_weights <- cr_info$sample_weights
    if (cr_info$has_censoring || cr_info$has_competing) {
      train_obs <- in_train & cr_info$observed_at_horizon
      if (sum(train_obs) < 20) {
        return(list(fold=fold, tau=rep(NA,sum(in_fold))))
      }
      Y_fg <- Y_event[train_obs]
      Y_fg[D[train_obs]==2] <- horizon+1
      forest <- causal_survival_forest(
        X[train_obs,,drop=FALSE], Y=Y_fg, D=D_event[train_obs], W=W[train_obs],
        sample.weights=sample_weights[train_obs], horizon=horizon
      )
    } else {
      Y_simple <- Y_event[in_train]
      Y_simple[D[in_train]==2] <- horizon+1
      forest <- causal_survival_forest(
        X[in_train,,drop=FALSE], Y=Y_simple, D=D_event[in_train], W=W[in_train],
        horizon=horizon
      )
    }
    tau_fold <- predict(forest, newdata=X[in_fold,,drop=FALSE])$predictions
    list(fold=fold, tau=tau_fold)
  }, .options=furrr_options(seed=TRUE))
  tau <- rep(NA,n_obs)
  for (res in tau_results) tau[fold_indices==res$fold] <- res$tau
  Y_pseudo <- Y
  Y_pseudo[D==2] <- horizon+1
  W_hat <- rep(0.5,n_obs)
  Y0_hat <- rep(mean(Y_pseudo[W==0]),n_obs)
  BLP_tests <- run_blp_tests_survival(Y_pseudo,W,W_hat,Y0_hat,tau,fold_indices,n_folds)
  BLP_whole <- BLP_survival(Y_pseudo,W,W_hat,Y0_hat,tau)$coefficients[,c(1,4)]
  indep_tests <- run_independence_tests(X,tau,fold_indices,n_folds)
  indep_whole <- run_independence_test_whole(X,tau)
  return(list(tau=data.frame(tau=tau),
              BLP_tests=BLP_tests,
              BLP_whole=BLP_whole,
              independence_tests=indep_tests,
              independence_whole=indep_whole,
              method="causal_survival_forest"))
}

###################
# Method 2: DR Oracle for Survival
###################
run_dr_oracle_survival <- function(X, Y, D, W, fmla_info, fold_indices, 
                                   fold_list, n_folds, horizon) {
  n_obs <- nrow(X)
  event1_params <- fmla_info$event1_params
  event2_params <- fmla_info$event2_params
  X_df <- as.data.frame(X)
  rmst_control <- numeric(n_obs)
  rmst_treat <- numeric(n_obs)
  for(i in 1:n_obs) {
    log_scale1 <- log(event1_params$scale_base) + 
      event1_params$b1_scale * X_df$X1[i] + 
      event1_params$b2_scale * X_df$X2[i]
    log_scale2 <- log(event2_params$scale_base) + 
      event2_params$b1_scale * X_df$X1[i] + 
      event2_params$b2_scale * X_df$X2[i]
    scale1 <- exp(log_scale1)
    scale2 <- exp(log_scale2)
    shape1_control <- event1_params$shape_base
    shape2_control <- event2_params$shape_base
    shape1_treat <- pmax(event1_params$shape_base,0.1)
    shape2_treat <- pmax(event2_params$shape_base,0.1)
    S_overall_control <- function(t) exp(-(t/scale1)^shape1_control-(t/scale2)^shape2_control)
    S_overall_treat <- function(t) exp(-(t/scale1)^shape1_treat-(t/scale2)^shape2_treat)
    lambda1_control <- function(t) (shape1_control/scale1)*(t/scale1)^(shape1_control-1)
    F1_control <- function(t) integrate(function(s)lambda1_control(s)*S_overall_control(s),0,t,rel.tol=1e-6)$value
    rmst_control[i] <- integrate(function(t)1-F1_control(t),0,horizon,rel.tol=1e-6)$value
    lambda1_treat <- function(t) (shape1_treat/scale1)*(t/scale1)^(shape1_treat-1)
    F1_treat <- function(t) integrate(function(s)lambda1_treat(s)*S_overall_treat(s),0,t,rel.tol=1e-6)$value
    rmst_treat[i] <- integrate(function(t)1-F1_treat(t),0,horizon,rel.tol=1e-6)$value
  }
  po <- rmst_treat - rmst_control
  tau <- estimate_cate_survival(X, po, fold_indices, fold_list, n_folds)
  Y_pseudo <- Y; Y_pseudo[D==2] <- horizon+1
  W_hat <- rep(0.5, n_obs)
  Y0_hat <- rmst_control
  BLP_tests <- run_blp_tests_survival(Y_pseudo,W,W_hat,Y0_hat,tau,fold_indices,n_folds)
  BLP_whole <- BLP_survival(Y_pseudo,W,W_hat,Y0_hat,tau)$coefficients[,c(1,4)]
  indep_tests <- run_independence_tests(X,tau,fold_indices,n_folds)
  indep_whole <- run_independence_test_whole(X,tau)
  return(list(
    tau = data.frame(tau = tau),
    BLP_tests = BLP_tests,
    BLP_whole = BLP_whole,
    independence_tests = indep_tests,
    independence_whole = indep_whole,
    method = "dr_oracle_survival"
  ))
}

###################
# Method 3: DR Semi-Oracle for Survival
###################
run_dr_semi_oracle_survival <- function(X, Y, D, W, fold_indices, 
                                        fold_list, n_folds, horizon) {
  n_obs <- nrow(X)
  W_hat <- rep(0.5, n_obs)
  cross_fits <- future_map(seq_len(n_folds), function(fold) {
    in_train <- !(fold_indices == fold)
    in_test <- which(!in_train)
    cr_info <- handle_censoring_and_competing_risks(X, Y, D, W, horizon, in_train)
    Y_event <- cr_info$Y_event
    D_event <- cr_info$D_event
    sample_weights <- cr_info$sample_weights
    observed_train <- in_train & cr_info$observed_at_horizon
    if (sum(observed_train) < 20) {
      return(list(po = rep(NA, length(in_test)), Y0_hat = rep(NA, length(in_test)), fold = fold))
    }
    # Estimate survival for event of interest under W=1 and W=0
    sf_event <- survival_forest(
      cbind(X[observed_train, , drop=FALSE], W[observed_train]),
      Y = Y_event[observed_train],
      D = as.integer(D[observed_train]==1),
      sample.weights = sample_weights[observed_train]
    )
    X_test0 <- cbind(X[in_test, , drop=FALSE], W=0)
    X_test1 <- cbind(X[in_test, , drop=FALSE], W=1)
    S1_0 <- predict(sf_event, newdata=X_test0, failure.times=horizon)$predictions
    S1_1 <- predict(sf_event, newdata=X_test1, failure.times=horizon)$predictions
    # RMST for event 1 under both treatments, by integration
    # (Numerically approximate RMST here for horizon)
    t_grid <- seq(0, horizon, length.out=101)
    get_rmst <- function(S_curve) sum(S_curve)*(horizon/100)
    S_curve0 <- predict(sf_event, newdata=X_test0, failure.times=t_grid)$predictions
    S_curve1 <- predict(sf_event, newdata=X_test1, failure.times=t_grid)$predictions
    rmst_0 <- apply(S_curve0, 1, get_rmst)
    rmst_1 <- apply(S_curve1, 1, get_rmst)
    cate <- rmst_1 - rmst_0
    # DR correction
    W_test <- W[in_test]
    observed_rmst <- as.numeric(D[in_test]==1 & Y[in_test]<=horizon)
    rmst_hat <- W_test*rmst_1 + (1-W_test)*rmst_0
    po <- cate + ((observed_rmst - rmst_hat)*(W_test-0.5))/(0.5*0.5)
    list(po=po, Y0_hat=rmst_0, fold=fold)
  }, .options=furrr_options(seed=TRUE))
  po <- rep(NA, n_obs)
  Y0_hat <- rep(NA, n_obs)
  for (i in seq_along(cross_fits)) {
    fold <- cross_fits[[i]]$fold
    in_fold <- fold_indices == fold
    po[in_fold] <- cross_fits[[i]]$po
    Y0_hat[in_fold] <- cross_fits[[i]]$Y0_hat
  }
  tau <- estimate_cate_survival(X, po, fold_indices, fold_list, n_folds)
  Y_pseudo <- Y; Y_pseudo[D==2] <- horizon+1
  BLP_tests <- run_blp_tests_survival(Y_pseudo,W,W_hat,Y0_hat,tau,fold_indices,n_folds)
  BLP_whole <- BLP_survival(Y_pseudo,W,W_hat,Y0_hat,tau)$coefficients[,c(1,4)]
  indep_tests <- run_independence_tests(X,tau,fold_indices,n_folds)
  indep_whole <- run_independence_test_whole(X,tau)
  return(list(
    tau = data.frame(tau = tau),
    BLP_tests = BLP_tests,
    BLP_whole = BLP_whole,
    independence_tests = indep_tests,
    independence_whole = indep_whole,
    method = "dr_semi_oracle_survival"
  ))
}

###################
# Method 4: Double-Robust Survival Forest (grf only)
###################
run_dr_survival_grf <- function(X, Y, D, W, fold_indices, fold_list, n_folds, horizon) {
  n_obs <- nrow(X)
  W_hat <- rep(0.5, n_obs)
  cross_fits <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_test <- which(!in_train)
    cr_info <- handle_censoring_and_competing_risks(X, Y, D, W, horizon, in_train)
    Y_event <- cr_info$Y_event
    D_event <- cr_info$D_event
    sample_weights <- cr_info$sample_weights
    observed_train <- in_train & cr_info$observed_at_horizon
    if (sum(observed_train) < 20) {
      return(list(po = rep(NA, length(in_test)), Y0_hat = rep(NA, length(in_test)), fold = fold))
    }
    sf_event <- survival_forest(
      cbind(X[observed_train, , drop=FALSE], W[observed_train]),
      Y = Y_event[observed_train],
      D = as.integer(D[observed_train]==1),
      sample.weights = sample_weights[observed_train]
    )
    X_test0 <- cbind(X[in_test, , drop=FALSE], W=0)
    X_test1 <- cbind(X[in_test, , drop=FALSE], W=1)
    S_curve0 <- predict(sf_event, newdata=X_test0, failure.times=seq(0, horizon, length.out=101))$predictions
    S_curve1 <- predict(sf_event, newdata=X_test1, failure.times=seq(0, horizon, length.out=101))$predictions
    get_rmst <- function(S_curve) sum(S_curve)*(horizon/100)
    rmst_0 <- apply(S_curve0, 1, get_rmst)
    rmst_1 <- apply(S_curve1, 1, get_rmst)
    cate <- rmst_1 - rmst_0
    # DR correction
    W_test <- W[in_test]
    observed_rmst <- as.numeric(D[in_test]==1 & Y[in_test]<=horizon)
    rmst_hat <- W_test*rmst_1 + (1-W_test)*rmst_0
    po <- cate + ((observed_rmst - rmst_hat)*(W_test-0.5))/(0.5*0.5)
    list(po=po, Y0_hat=rmst_0, fold=fold)
  }, .options=furrr_options(seed=TRUE))
  po <- rep(NA, n_obs)
  Y0_hat <- rep(NA, n_obs)
  for (i in seq_along(cross_fits)) {
    fold <- cross_fits[[i]]$fold
    in_fold <- fold_indices == fold
    po[in_fold] <- cross_fits[[i]]$po
    Y0_hat[in_fold] <- cross_fits[[i]]$Y0_hat
  }
  tau <- estimate_cate_survival(X, po, fold_indices, fold_list, n_folds)
  Y_pseudo <- Y; Y_pseudo[D==2] <- horizon+1
  BLP_tests <- run_blp_tests_survival(Y_pseudo,W,W_hat,Y0_hat,tau,fold_indices,n_folds)
  BLP_whole <- BLP_survival(Y_pseudo,W,W_hat,Y0_hat,tau)$coefficients[,c(1,4)]
  indep_tests <- run_independence_tests(X,tau,fold_indices,n_folds)
  indep_whole <- run_independence_test_whole(X,tau)
  return(list(
    tau = data.frame(tau = tau),
    BLP_tests = BLP_tests,
    BLP_whole = BLP_whole,
    independence_tests = indep_tests,
    independence_whole = indep_whole,
    method = "dr_survival_grf"
  ))
}

###################
# Helper Functions for CATE estimation and tests
###################
estimate_cate_survival <- function(X, po, fold_indices, fold_list, n_folds) {
  n_obs <- nrow(X)
  tau_results <- future_map(seq_along(fold_list), function(i) {
    fold <- fold_list[i]
    in_train <- fold_indices != fold
    in_fold <- !in_train
    if (sum(!is.na(po[in_train])) < 10) {
      return(list(fold = fold, predictions = rep(NA, sum(in_fold))))
    }
    forest <- regression_forest(X[in_train, , drop = FALSE], po[in_train])
    tau_pred <- predict(forest, newdata = X[in_fold, , drop = FALSE])$predictions
    list(fold = fold, predictions = tau_pred)
  }, .options = furrr_options(seed = TRUE))
  tau <- rep(NA, n_obs)
  for (result in tau_results) {
    fold <- result$fold
    in_fold <- fold_indices == fold
    tau[in_fold] <- result$predictions
  }
  return(tau)
}

BLP_survival <- function(Y,W,W_hat,Y0_hat,tau){
  cc <- complete.cases(Y,W,W_hat,Y0_hat,tau)
  if(sum(cc)<10) return(list(coefficients=matrix(NA,nrow=4,ncol=4)))
  Y <- Y[cc]; W <- W[cc]; W_hat <- W_hat[cc]; Y0_hat <- Y0_hat[cc]; tau <- tau[cc]
  tau_c <- tau - mean(tau,na.rm=TRUE)
  df <- data.frame(Y=Y,
                   W=W,
                   Wdiff=W-W_hat,
                   tauTerm=tau_c*(W-W_hat))
  lm_fit <- lm(Y ~ W + Wdiff + tauTerm, data=df)
  return(summary(lm_fit))
}

run_blp_tests_survival <- function(Y, W, W_hat, Y0_hat, tau, fold_indices, n_folds) {
  future_map(seq_len(n_folds), function(fold){
    in_test <- fold_indices==fold
    BLP_res <- BLP_survival(Y[in_test],W[in_test],W_hat[in_test],Y0_hat[in_test],tau[in_test])
    if(!is.null(BLP_res$coefficients) && nrow(BLP_res$coefficients)>=3){
      data.frame(fold=fold,
                 beta1_est=BLP_res$coefficients[2,1],
                 beta1_pval=BLP_res$coefficients[2,4],
                 beta2_est=BLP_res$coefficients[3,1],
                 beta2_pval=BLP_res$coefficients[3,4])
    } else data.frame(fold=fold,beta1_est=NA,beta1_pval=NA,beta2_est=NA,beta2_pval=NA)
  },.options=furrr_options(seed=TRUE)) %>% bind_rows()
}

run_independence_tests <- function(X,tau,fold_indices,n_folds){
  future_map(seq_len(n_folds), function(fold){
    in_test <- fold_indices==fold
    if(sum(!is.na(tau[in_test]))<10) return(data.frame(fold=fold,pval=NA))
    df <- data.frame(tau=tau[in_test], X[in_test,,drop=FALSE])
    itest <- coin::independence_test(tau ~ ., data=df)
    data.frame(fold=fold,pval=coin::pvalue(itest))
  },.options=furrr_options(seed=TRUE)) %>% bind_rows()
}
run_independence_test_whole <- function(X,tau){
  if(sum(!is.na(tau))<10) return(list(pvalue=NA))
  df <- data.frame(tau=tau[!is.na(tau)], X[!is.na(tau),,drop=FALSE])
  itest <- coin::independence_test(tau ~ ., data=df)
  list(pvalue=coin::pvalue(itest))
}
