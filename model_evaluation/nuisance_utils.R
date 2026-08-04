##########
# nuisance estimation functions
##########
# general
require(h2o)
require(xgboost)
unlist_order <- function(reslist) {
  df <- do.call(rbind, reslist)
  df <- df %>%
    arrange(row_index) %>%
    select(-row_index)
  return(df)
}

calculate_pseudos <- function(df, Y, W) {
  df <- df %>%
    mutate(
      mu_DR = mu1_DR * W + (1 - W) * mu0_DR
    )
  df <- df %>%
    mutate(
      tau_T = mu1_T - mu0_T,
      phi = mu1_DR - mu0_DR + ((Y - mu_DR) * (W - pi)) / (pi * (1 - pi)),
      phi05 = mu1_DR - mu0_DR + ((Y - mu_DR) * (W - 0.5)) / (0.5 * (1 - 0.5))
    )
  return(df)
}

# XGB
create_param_grid <- function(n_cores = 1) {
  grid_base <- expand.grid(
    nthread = n_cores,
    eta = c(0.1, 0.2, 0.3),
    max_depth = c(2, 4, 6),
    subsample = c(0.6, 0.8),
    colsample_bytree = c(0.7, 1)
  )

  list(
    continuous = grid_base %>% mutate(objective = "reg:squarederror"),
    binary = grid_base %>% mutate(objective = "binary:logistic")
  )
}

prepare_xgb_matrices <- function(X, Y, W, filter_indices, training = TRUE) {
  if (training) {
    list(
      X_Y = xgb.DMatrix(
        data = as.matrix(X[filter_indices, ]),
        label = Y[filter_indices]
      ),
      X_0_Y = xgb.DMatrix(
        data = as.matrix(X[filter_indices & W == 0, ]),
        label = Y[filter_indices & W == 0]
      ),
      X_1_Y = xgb.DMatrix(
        data = as.matrix(X[filter_indices & W == 1, ]),
        label = Y[filter_indices & W == 1]
      ),
      XW_Y = xgb.DMatrix(
        data = as.matrix(cbind(W = W[filter_indices], X[filter_indices, ])),
        label = Y[filter_indices]
      ),
      X_W = xgb.DMatrix(
        data = as.matrix(X[filter_indices, ]),
        label = W[filter_indices]
      )
    )
  } else {
    # For prediction matrices
    list(
      X = xgb.DMatrix(data = as.matrix(X[filter_indices, ])),
      X0 = xgb.DMatrix(data = as.matrix(cbind(W = 0, X[filter_indices, ]))),
      X1 = xgb.DMatrix(data = as.matrix(cbind(W = 1, X[filter_indices, ])))
    )
  }
}

fit_xgb <- function(grid, data_matrix) {
  cv_results <- run_xgb_cv(grid, data_matrix)
  trained_result <- train_best_xgb(grid, cv_results, data_matrix)
  return(trained_result$model)
}

run_xgb_cv <- function(grid, train_data) {
  cv_results <- apply(grid, 1, function(params) {
    xgb.cv(
      params = as.list(params),
      data = train_data,
      nrounds = 100,
      nfold = 5,
      early_stopping_rounds = 10,
      verbose = 0
    )
  })

  cv_results <- lapply(cv_results, function(params) {
    error_col <- if (grid$objective[1] == "reg:squarederror") {
      "test_rmse_mean"
    } else {
      "test_logloss_mean"
    }
    data.frame(
      nrounds = params$best_iteration,
      error = min(params$evaluation_log[[error_col]])
    )
  })

  return(do.call(rbind, cv_results))
}

train_best_xgb <- function(grid, cv_results, train_data) {
  best_idx <- which.min(cv_results$error)
  best_params <- as.list(grid[best_idx, ])
  best_nrounds <- cv_results$nrounds[best_idx]

  best_model <- xgb.train(
    params = best_params,
    data = train_data,
    nrounds = best_nrounds,
    verbose = 0
  )

  return(list(model = best_model, params = best_params, nrounds = best_nrounds))
}

fit_nuisance_xgb <- function(train_matrices, param_grids) {
  list(
    mu0_T = fit_xgb(param_grids$continuous, train_matrices$X_0_Y),
    mu1_T = fit_xgb(param_grids$continuous, train_matrices$X_1_Y),
    mu_DR = fit_xgb(param_grids$continuous, train_matrices$XW_Y),
    pi = fit_xgb(param_grids$binary, train_matrices$X_W)
  )
}

predict_nuisance_xgb <- function(models, pred_matrices) {
  data.frame(
    mu0_T = predict(models$mu0_T, newdata = pred_matrices$X),
    mu1_T = predict(models$mu1_T, newdata = pred_matrices$X),
    mu0_DR = predict(models$mu_DR, newdata = pred_matrices$X0),
    mu1_DR = predict(models$mu_DR, newdata = pred_matrices$X1),
    pi = predict(models$pi, newdata = pred_matrices$X)
  )
}

run_xgb_cross_validation <- function(
  X,
  Y,
  W,
  fold_indices,
  fold_list,
  param_grids,
  method = "cv"
) {
  nuisances <- lapply(fold_list, function(fold) {
    # Determine train/test filters based on method
    if (method == "cv") {
      train_filter <- !(fold_indices %in% fold)
      test_filter <- !train_filter
    } else if (method == "infold") {
      train_filter <- fold_indices == fold
      test_filter <- train_filter
    }

    # Prepare data matrices
    train_matrices <- prepare_xgb_matrices(
      X,
      Y,
      W,
      train_filter,
      training = TRUE
    )
    pred_matrices <- prepare_xgb_matrices(
      X,
      Y,
      W,
      test_filter,
      training = FALSE
    )

    # Fit models
    models <- fit_nuisance_xgb(train_matrices, param_grids)

    # Make predictions
    predictions <- predict_nuisance_xgb(models, pred_matrices)
    predictions$row_index <- which(test_filter)

    return(predictions)
  })

  # Combine results and calculate pseudo-outcomes
  result <- unlist_order(nuisances)
  result <- calculate_pseudos(result, Y, W)

  return(result)
}

run_xgb_whole_dataset <- function(X, Y, W, param_grids) {
  n_obs <- nrow(X)
  all_indices <- rep(TRUE, n_obs)

  # Prepare data matrices
  train_matrices <- prepare_xgb_matrices(X, Y, W, all_indices, training = TRUE)
  pred_matrices <- prepare_xgb_matrices(X, Y, W, all_indices, training = FALSE)

  # Fit models
  models <- fit_nuisance_xgb(train_matrices, param_grids)

  # Make predictions
  result <- predict_nuisance_xgb(models, pred_matrices)
  result <- calculate_pseudos(result, Y, W)

  return(result)
}

run_all_xgb_nuisance <- function(
  X,
  Y,
  W,
  fold_indices,
  fold_list,
  n_cores = 1
) {
  param_grids <- create_param_grid(n_cores)

  # Cross-validation
  xgb_cv <- run_xgb_cross_validation(
    X,
    Y,
    W,
    fold_indices,
    fold_list,
    param_grids,
    "cv"
  )

  # In-fold
  xgb_infold <- run_xgb_cross_validation(
    X,
    Y,
    W,
    fold_indices,
    fold_list,
    param_grids,
    "infold"
  )

  # Whole dataset
  xgb_whole <- run_xgb_whole_dataset(X, Y, W, param_grids)

  return(list(cv = xgb_cv, infold = xgb_infold, whole = xgb_whole))
}

# AutoML
setup_h2o_cluster <- function(n_cores, mem, model_seed = NULL) {
  # Check if cluster is already running
  tryCatch(
    {
      h2o.clusterIsUp()
      cat("H2O cluster is already running\n")
    },
    error = function(e) {
      # Initialize H2O
      h2o.init(
        port = parallelly::freePort(),
        nthreads = n_cores,
        max_mem_size = mem,
        min_mem_size = "1G"
      )
    }
  )
}

# Set up shutdown handler
h2o_shutdown_check <- function() {
  cat("Ensuring H2O cluster shutdown...\n")
  tryCatch(
    {
      h2o.shutdown(prompt = FALSE)
      cat("H2O cluster shutdown successfully\n")
    },
    error = function(e) {
      cat("H2O cluster was already shutdown\n")
    }
  )
}

prepare_h2o_data <- function(X, Y, W, filter_indices, training = TRUE) {
  # Prepare different data combinations
  if (training) {
    YX <- cbind(Y, X)
    YWX <- cbind(Y, W, X)
    WX <- cbind(W = as.factor(W), X)

    list(
      YX_0 = as.h2o(YX[filter_indices & W == 0, ]),
      YX_1 = as.h2o(YX[filter_indices & W == 1, ]),
      YWX = as.h2o(YWX[filter_indices, ]),
      WX = as.h2o(WX[filter_indices, ])
    )
  } else {
    list(
      X = as.h2o(X[filter_indices, ]),
      X0 = as.h2o(cbind(W = 0, X)[filter_indices, ]),
      X1 = as.h2o(cbind(W = 1, X)[filter_indices, ])
    )
  }
}

fit_automl_models <- function(
  train_data,
  model_seed,
  max_models = 20,
  exclude_algos = c("DeepLearning", "XGBoost")
) {
  list(
    mu0_T = h2o.automl(
      y = "Y",
      training_frame = train_data$YX_0,
      nfolds = 0,
      max_models = max_models,
      exclude_algos = exclude_algos,
      seed = model_seed
    ),
    mu1_T = h2o.automl(
      y = "Y",
      training_frame = train_data$YX_1,
      nfolds = 0,
      max_models = max_models,
      exclude_algos = exclude_algos,
      seed = model_seed
    ),
    mu_DR = h2o.automl(
      y = "Y",
      training_frame = train_data$YWX,
      nfolds = 0,
      max_models = max_models,
      exclude_algos = exclude_algos,
      seed = model_seed
    ),
    pi = h2o.automl(
      y = "W",
      training_frame = train_data$WX,
      nfolds = 0,
      max_models = max_models,
      exclude_algos = exclude_algos,
      seed = model_seed
    )
  )
}

predict_automl_models <- function(models, pred_data) {
  pi_pred <- h2o.predict(models$pi, newdata = pred_data$X) %>% as.data.frame()

  data.frame(
    mu0_T = h2o.predict(models$mu0_T, newdata = pred_data$X) %>% as.vector(),
    mu1_T = h2o.predict(models$mu1_T, newdata = pred_data$X) %>% as.vector(),
    mu0_DR = h2o.predict(models$mu_DR, newdata = pred_data$X0) %>% as.vector(),
    mu1_DR = h2o.predict(models$mu_DR, newdata = pred_data$X1) %>% as.vector(),
    pi = pi_pred$p1
  )
}

run_automl_cross_validation <- function(
  X,
  Y,
  W,
  fold_indices,
  fold_list,
  model_seed,
  method = "cv",
  max_models = 20,
  exclude_algos = c("DeepLearning", "XGBoost")
) {
  nuisances <- lapply(fold_list, function(fold) {
    # Determine train/test filters based on method
    if (method == "cv") {
      train_filter <- !(fold_indices %in% fold)
      test_filter <- !train_filter
    } else if (method == "infold") {
      train_filter <- fold_indices == fold
      test_filter <- train_filter
    }

    # Prepare data
    train_data <- prepare_h2o_data(X, Y, W, train_filter, training = TRUE)
    pred_data <- prepare_h2o_data(X, Y, W, test_filter, training = FALSE)

    # Fit models
    models <- fit_automl_models(
      train_data,
      model_seed,
      max_models,
      exclude_algos
    )

    # Make predictions
    predictions <- predict_automl_models(models, pred_data)
    predictions$row_index <- which(test_filter)

    # Clean up H2O objects for this fold
    h2o.removeAll()

    return(predictions)
  })

  # Combine results and calculate pseudo-outcomes
  result <- unlist_order(nuisances)

  # Handle W conversion for infold method (based on your original code)
  W_final <- if (method == "infold") as.numeric(W) - 1 else W
  result <- calculate_pseudos(result, Y, W_final)

  return(result)
}

run_automl_whole_dataset <- function(
  X,
  Y,
  W,
  model_seed,
  max_models = 20,
  exclude_algos = c("DeepLearning", "XGBoost")
) {
  n_obs <- nrow(X)
  all_indices <- rep(TRUE, n_obs)

  train_data <- prepare_h2o_data(X, Y, W, all_indices, training = TRUE)
  pred_data <- prepare_h2o_data(X, Y, W, all_indices, training = FALSE)

  # Fit models
  models <- fit_automl_models(train_data, model_seed, max_models, exclude_algos)

  # Prepare prediction data

  # Make predictions
  result <- predict_automl_models(models, pred_data)
  result <- calculate_pseudos(result, Y, W)

  return(result)
}

run_all_automl_nuisance <- function(
  X,
  Y,
  W,
  fold_indices,
  fold_list,
  n_cores = 1,
  mem = "5G",
  model_seed,
  max_models = 20,
  exclude_algos = c("DeepLearning", "XGBoost")
) {
  # Setup H2O cluster
  setup_h2o_cluster(n_cores = n_cores, mem = mem, model_seed = model_seed)
  on.exit(h2o_shutdown_check, add = TRUE)

  # Cross-validation
  auto_cv <- run_automl_cross_validation(
    X,
    Y,
    W,
    fold_indices,
    fold_list,
    model_seed,
    "cv",
    max_models,
    exclude_algos
  )

  # In-fold
  auto_infold <- run_automl_cross_validation(
    X,
    Y,
    W,
    fold_indices,
    fold_list,
    model_seed,
    "infold",
    max_models,
    exclude_algos
  )

  # Whole dataset
  auto_whole <- run_automl_whole_dataset(
    X,
    Y,
    W,
    model_seed,
    max_models,
    exclude_algos
  )

  # Shutdown H2O
  h2o.shutdown(prompt = FALSE)

  return(list(cv = auto_cv, infold = auto_infold, whole = auto_whole))
}
