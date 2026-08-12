###############
# title: data-processing helpers - model evaluation study
###############
# setup_rng_seed() (the old local RNG-stream reimplementation) and
# generate_data() (the old benchtm call) are gone: setup_rng_stream() now
# comes from R/utils.R (statement-for-statement identical to the old
# setup_rng_seed(), just named consistently with every other study), and
# data generation comes from me_dgms.R.
#
# prepare_design_matrix()/split_folds() are otherwise unchanged from the
# original sim_utils.R. Both select columns by name, not position, so they
# were already robust to the DGM swap - confirmed empirically, not just by
# inspection: running R/dgm_scenarios.R's output through this exact
# model.matrix()/preProcess() logic reproduces the covariate columns exactly,
# since every covariate the new DGM emits is already numeric (no factors
# left to expand).

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
#
# No longer returns fold_pairs: me_models.R moved off double crossfitting
# (fold pairs, C(V,2) fits) onto single crossfitting (fold_list, V fits) - see
# me_models.R's header. fold_pairs was only ever consumed there.
split_folds <- function(Y, k = 10) {
  fold_indices <- caret::createFolds(y = Y, k = k, list = F)
  fold_list <- sort(unique(fold_indices))
  list(
    fold_indices = fold_indices,
    fold_list = fold_list
  )
}
