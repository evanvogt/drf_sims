##########
# model fitting functions
##########
require(future.apply)
require(ranger)
require(SuperLearner)
# hyper parameter specifications for CATE learners
create_rf_hyperparams <- function(
  mtry = NULL,
  max.depth = NUL,
  sample.fraction = 0.5,
  replace = FALSE
) {
  list(
    mtry = mtry,
    max.depth = max.depth,
    sample.fraction = sample.fraction,
    replace = replace
  )
}

create_net_hyperparams <- function(alpha) {
  list(alpha = alpha)
}

create_SL_hyperparams <- function(method, SL.library) {
  list(method = method, SL.library = SL.library)
}

# general fitting functions
collate_predictions <- function(
  fold_list,
  fold_pairs,
  fold_indices,
  reslist,
  target
) {
  lapply(fold_list, function(fold) {
    predictions <- rep(NA, length(fold_indices))
    for (j in seq_along(fold_pairs)) {
      if (fold %in% fold_pairs[[j]]) {
        predictions[fold_indices %in% fold_pairs[[j]]] <- reslist[[j]][[target]]
      }
    }
    predictions[fold_indices == fold] <- NA
    predictions
  }) %>%
    simplify2array()
}

collect_tau <- function(tau_list, fold_list, fold_indices) {
  tau <- rep(NA, length(fold_indices))
  for (fold in fold_list) {
    tau[fold_indices == fold] <- tau_list[[fold]]
  }
  return(tau)
}
# fitting functions for the random forest
fit_rf <- function(
  Y,
  W,
  X,
  hyper_list = list(),
  fold_indices,
  fold_list,
  fold_pairs
) {
  hyper_list$write.forest <- TRUE

  nuisances <- future_lapply(
    fold_pairs,
    crossfit_double_RF,
    Y,
    W,
    X,
    hyper_list,
    fold_indices,
    future.seed = TRUE
  )

  predictions <- c("po", "Y0.hat", "Y1.hat", "W.hat")
  matrices <- setNames(
    lapply(predictions, function(pred) {
      collate_predictions(fold_list, fold_pairs, fold_indices, nuisances, pred)
    }),
    predictions
  )

  tau <- fit_tau_rf(X, fold_list, fold_indices, matrices$po, hyper_list)

  c(list(tau = tau), matrices)
}

crossfit_double_RF <- function(fold_pair, Y, W, X, hyper_list, fold_indices) {
  train_filter <- !(fold_indices %in% fold_pair)
  test_filter <- !train_filter

  fit_model <- function(formula, data) {
    hyper <- hyper_list
    hyper$data <- data
    hyper$formula <- formula
    do.call(ranger, hyper)
  }

  Y_model <- fit_model(
    Y ~ .,
    data.frame(Y = Y[train_filter], W = W[train_filter], X[train_filter, ])
  )
  W_model <- fit_model(
    W ~ .,
    data.frame(W = W[train_filter], X[train_filter, ])
  )

  test_data <- data.frame(W, X)[test_filter, ]
  Y0.hat <- predict(Y_model, mutate(test_data, W = 0))$predictions
  Y1.hat <- predict(Y_model, mutate(test_data, W = 1))$predictions
  W.hat <- predict(W_model, X[test_filter, ])$predictions

  W_test <- W[test_filter]
  Y_test <- Y[test_filter]
  Y.hat <- W_test * Y1.hat + (1 - W_test) * Y0.hat

  po <- Y1.hat -
    Y0.hat +
    ((Y_test - Y.hat) * (W_test - W.hat)) / (W.hat * (1 - W.hat))

  list(po = po, Y0.hat = Y0.hat, Y1.hat = Y1.hat, W.hat = W.hat)
}

fit_tau_rf <- function(X, fold_list, fold_indices, po_matrix, hyper_list) {
  tau <- rep(NA, length(fold_indices))

  tau_list <- future_lapply(
    fold_list,
    function(fold) {
      train_filter <- !(fold_indices %in% fold)
      test_filter <- !train_filter

      hyper_po <- hyper_list
      hyper_po$data <- data.frame(
        po = po_matrix[train_filter, fold],
        X[train_filter, ]
      )
      hyper_po$formula <- po ~ .

      po_model <- do.call(ranger, hyper_po)
      tau[test_filter] <- predict(po_model, X[test_filter, ])$predictions
    },
    future.seed = TRUE
  )

  for (fold in fold_list) {
    tau[fold_indices == fold] <- tau_list[[fold]]
  }

  tau
}

# fitting functions for GLMnet
fit_glmnet <- function(
  Y,
  W,
  X,
  hyper_list = list(),
  fold_indices,
  fold_list,
  fold_pairs
) {
  nuisances <- future_lapply(
    fold_pairs,
    crossfit_double_glmnet,
    Y,
    W,
    X,
    hyper_list,
    fold_indices,
    future.seed = TRUE
  )

  predictions <- c("po", "Y0.hat", "Y1.hat", "W.hat")
  matrices <- setNames(
    lapply(predictions, function(pred) {
      collate_predictions(fold_list, fold_pairs, fold_indices, nuisances, pred)
    }),
    predictions
  )

  tau <- fit_tau_glmnet(X, fold_list, fold_indices, matrices$po, hyper_list)

  c(list(tau = tau), matrices)
}

crossfit_double_glmnet <- function(
  fold_pair,
  Y,
  W,
  X,
  hyper_list,
  fold_indices
) {
  custom_net <- create.Learner("SL.glmnet", params = hyper_list)

  train_filter <- !(fold_indices %in% fold_pair)
  test_filter <- !train_filter

  fit_model <- function(y, x, family) {
    SuperLearner(
      y,
      x,
      family = family,
      SL.library = custom_net$names,
      method = "method.CC_LS"
    )
  }

  Y_model <- fit_model(Y[train_filter], cbind(W, X)[train_filter, ], gaussian())
  W_model <- fit_model(W[train_filter], X[train_filter, ], binomial())

  test_data <- data.frame(W, X)[test_filter, ]
  Y0.hat <- predict(Y_model, mutate(test_data, W = 0))$pred
  Y1.hat <- predict(Y_model, mutate(test_data, W = 1))$pred
  W.hat <- predict(W_model, X[test_filter, ])$pred

  W_test <- W[test_filter]
  Y_test <- Y[test_filter]
  Y.hat <- W_test * Y1.hat + (1 - W_test) * Y0.hat

  po <- Y1.hat -
    Y0.hat +
    ((Y_test - Y.hat) * (W_test - W.hat)) / (W.hat * (1 - W.hat))

  list(po = po, Y0.hat = Y0.hat, Y1.hat = Y1.hat, W.hat = W.hat)
}

fit_tau_glmnet <- function(X, fold_list, fold_indices, po_matrix, hyper_list) {
  custom_net <- create.Learner("SL.glmnet", params = hyper_list)

  tau_list <- future_lapply(
    fold_list,
    function(fold) {
      train_filter <- !(fold_indices %in% fold)
      test_filter <- !train_filter

      po_model <- SuperLearner(
        po_matrix[train_filter, fold],
        X[train_filter, ],
        family = gaussian(),
        SL.library = custom_net$names,
        method = "method.CC_LS"
      )

      predict(po_model, X[test_filter, ])$pred
    },
    future.seed = TRUE
  )

  collect_tau(tau_list, fold_list, fold_indices)
}

# fitting functions for SuperLearner
fit_SL <- function(
  Y,
  W,
  X,
  hyper_list = list(),
  fold_indices,
  fold_list,
  fold_pairs
) {
  nuisances <- future_lapply(
    fold_pairs,
    crossfit_double_SL,
    Y,
    W,
    X,
    hyper_list,
    fold_indices,
    future.seed = TRUE
  )

  predictions <- c("po", "Y0.hat", "Y1.hat", "W.hat")
  matrices <- setNames(
    lapply(predictions, function(pred) {
      collate_predictions(fold_list, fold_pairs, fold_indices, nuisances, pred)
    }),
    predictions
  )

  tau <- fit_tau_SL(X, fold_list, fold_indices, matrices$po, hyper_list)

  c(list(tau = tau), matrices)
}

crossfit_double_SL <- function(fold_pair, Y, W, X, hyper_list, fold_indices) {
  train_filter <- !(fold_indices %in% fold_pair)
  test_filter <- !train_filter

  fit_model <- function(y, x, family) {
    hyper <- hyper_list
    hyper$family <- family
    hyper$Y <- y
    hyper$X <- x
    do.call(SuperLearner, hyper)
  }

  Y_model <- fit_model(Y[train_filter], cbind(W, X)[train_filter, ], gaussian())
  W_model <- fit_model(W[train_filter], X[train_filter, ], binomial())

  test_data <- data.frame(W, X)[test_filter, ]
  Y0.hat <- predict(Y_model, mutate(test_data, W = 0))$pred
  Y1.hat <- predict(Y_model, mutate(test_data, W = 1))$pred
  W.hat <- predict(W_model, X[test_filter, ])$pred

  W_test <- W[test_filter]
  Y_test <- Y[test_filter]
  Y.hat <- W_test * Y1.hat + (1 - W_test) * Y0.hat

  po <- Y1.hat -
    Y0.hat +
    ((Y_test - Y.hat) * (W_test - W.hat)) / (W.hat * (1 - W.hat))

  list(po = po, Y0.hat = Y0.hat, Y1.hat = Y1.hat, W.hat = W.hat)
}

fit_tau_SL <- function(X, fold_list, fold_indices, po_matrix, hyper_list) {
  tau_list <- future_lapply(
    fold_list,
    function(fold) {
      train_filter <- !(fold_indices %in% fold)
      test_filter <- !train_filter

      hyper_po <- hyper_list
      hyper_po$family <- gaussian()
      hyper_po$X <- X[train_filter, ]
      hyper_po$Y <- po_matrix[train_filter, fold]

      po_model <- do.call(SuperLearner, hyper_po)
      predict(po_model, X[test_filter, ])$pred
    },
    future.seed = TRUE
  )

  tau <- rep(NA, length(fold_indices))
  for (fold in fold_list) {
    tau[fold_indices == fold] <- tau_list[[fold]]
  }

  tau
}
