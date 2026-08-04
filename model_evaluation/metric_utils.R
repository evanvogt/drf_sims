##########
# metric calculation functions
##########

# Extract surrogate values from nuisance datafiles
extract_surrogate <- function(nuisance, var_name) {
  sapply(nuisance, function(df) df[[var_name]], simplify = TRUE)
}

# Function to calculate scores for a given configuration and fold type
calc_scores_for_config <- function(
  config,
  fold_type,
  candidate_preds,
  surrogates,
  Y,
  W
) {
  if (config$type == "infl") {
    tau_surr <- surrogates[[config$surr_tau]][, fold_type]
    pi_surr <- surrogates[[config$surr_pi]][, fold_type]
    apply(candidate_preds, 2, function(tau_pred) {
      calc_infl_score(tau_pred, tau_surr, pi_surr, Y, W)
    })
  } else {
    # "dr" case
    phi_surr <- surrogates[[config$surr_phi]][, fold_type]
    apply(candidate_preds, 2, function(tau_pred) {
      calc_dr_risk(tau_pred, phi_surr)
    })
  }
}

# Influence corrected loss function
calc_infl_score <- function(tau_hat, tau, pi, Y, W) {
  A <- W - pi
  C <- pi * (1 - pi)
  B <- 2 * W * A / C
  diff_tau <- tau_hat - tau

  pehe_infl <- (1 - B) * tau^2 + B * Y * diff_tau - A * diff_tau^2 + tau_hat^2

  mean(pehe_infl)
}

# DR risk
calc_dr_risk <- function(tau_hat, tau_dr) {
  mean((tau_hat - tau_dr)^2)
}

# Fix automl naming issue
fix_automl <- function(nuisance_models) {
  for (fold_type in names(nuisance_models[["automl"]])) {
    names_vec <- names(nuisance_models[["automl"]][[fold_type]])
    if (length(names_vec) >= 5 && names_vec[5] == "p1") {
      names_vec[5] <- "pi"
      names(nuisance_models[["automl"]][[fold_type]]) <- names_vec
    }
  }
  nuisance_models
}

# Main metrics calculation function
calc_metrics <- function(
  models_list,
  nuisance_models,
  Y,
  W,
  truth_avail = FALSE,
  true_tau = NULL
) {
  # Added Y, W parameters
  # fix automl naming
  nuisance_models <- fix_automl(nuisance_models)

  # Format model predictions
  candidate_preds <- sapply(models_list, function(fit) fit$tau, simplify = TRUE)

  # Extract surrogates
  model_names <- names(nuisance_models)
  surrogates <- list()

  for (model_name in model_names) {
    model_data <- nuisance_models[[model_name]]
    surrogates[[paste0(
      "tau_T_",
      model_name
    )]] <- as.data.frame(extract_surrogate(model_data, "tau_T"))
    surrogates[[paste0("pi_", model_name)]] <- as.data.frame(extract_surrogate(
      model_data,
      "pi"
    ))
    surrogates[[paste0("phi_", model_name)]] <- as.data.frame(extract_surrogate(
      model_data,
      "phi"
    ))
  }

  # Generate score configurations
  score_configs <- unlist(
    lapply(model_names, function(model_name) {
      list(
        list(
          type = "infl",
          surr_tau = paste0("tau_T_", model_name),
          surr_pi = paste0("pi_", model_name),
          suffix = paste0("_", model_name)
        ),
        list(
          type = "dr",
          surr_phi = paste0("phi_", model_name),
          suffix = paste0("_", model_name)
        )
      )
    }),
    recursive = FALSE
  )

  fold_types <- names(nuisance_models[[model_names[1]]]) # Use first model for fold types

  # Initialize results dataframe
  metrics_all <- data.frame(model = colnames(candidate_preds))

  # Calculate scores for all configurations and fold types
  for (fold_type in fold_types) {
    for (config in score_configs) {
      scores <- calc_scores_for_config(
        config,
        fold_type,
        candidate_preds,
        surrogates,
        Y,
        W
      )
      col_name <- paste0(config$type, "_", fold_type, config$suffix)
      metrics_all[[col_name]] <- scores
    }
  }

  # Calculate true PEHE
  if (truth_avail) {
    metrics_all[["true_pehe"]] <- apply(candidate_preds, 2, function(tau_pred) {
      mean((tau_pred - true_tau)^2)
    })
  }

  metrics_all
}
