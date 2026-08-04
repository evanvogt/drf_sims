###############
# Simulation specific functions
###############

# uncorrelated random number generation
setup_rng_seed <- function(sim_num, base_seed = 1998) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(base_seed)

  # new stream for each simulation run
  for (i in 1:sim_num) {
    .Random.seed <<- parallel::nextRNGStream(.Random.seed)
  }
}

# Function to generate data
generate_data <- function(scen_param, scen_index) {
  scens <- scen_param %>% filter(type == "continuous")
  adat <- generate_scen_data(scen = scens[scen_index, ], include_truth = TRUE)
  true_tau <- adat$trt_effect
  adat <- adat %>%
    rename(W = trt) %>%
    select(-c(trt_effect, prob_diff))
  list(adat = adat, true_tau = true_tau)
}

# Function to process data
prepare_design_matrix <- function(data) {
  Y <- data$Y
  W <- data$W
  X <- data %>% select(-c(Y, W))
  covariates <- names(X)
  fmla <- formula(paste0("~ ", paste0(covariates, collapse = " + ")))
  X <- model.matrix(fmla, X)
  X <- data.frame(X[, -1]) # Remove intercept column
  scaling <- caret::preProcess(X)
  X <- predict(scaling, X)
  list(Y = Y, W = W, X = X, covariates = covariates, scaling = scaling)
}

# Function to split folds
split_folds <- function(Y, k = 10) {
  fold_indices <- caret::createFolds(y = Y, k = k, list = F)
  fold_list <- sort(unique(fold_indices))
  fold_pairs <- utils::combn(fold_list, 2, simplify = FALSE)
  list(
    fold_indices = fold_indices,
    fold_list = fold_list,
    fold_pairs = fold_pairs
  )
}
