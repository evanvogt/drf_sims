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

#' Blocks for the `holdout` nuisance arm, derived from the candidate's folds
#'
#' The `holdout` arm fits its nuisance on nothing but the fold the candidate
#' held out, so its blocks MUST come from the candidate's own fold assignment -
#' a fresh draw of the same size would not be "the data the candidate never
#' saw". This derives them rather than drawing them, which is also why it needs
#' no RNG and can run from a saved fold_info.
#'
#' Where a fold is too small to fit the two nuisance targets in (mu_DR and pi,
#' both fit on the whole block), adjacent folds are POOLED. At n = 250,
#' V = 10 that is 25-row folds pooled into 5 blocks of 50. The cost of
#' pooling, worth stating: a pooled block is only half-decoupled from any
#' single candidate fold model, since the partner fold WAS in that model's
#' training set. The alternative - a 25-row block - is not estimable at all,
#' so this is the lesser problem.
#'
#' @param fold_indices the candidate fold assignment (fold_info$fold_indices)
#' @param min_block smallest acceptable block size, default HOLDOUT_MIN_BLOCK
#' @return the same shape as split_folds(): fold_indices and fold_list
holdout_blocks <- function(fold_indices, min_block = HOLDOUT_MIN_BLOCK) {
  n_obs <- length(fold_indices)
  n_folds <- length(unique(fold_indices))
  block_size <- n_obs / n_folds

  pool <- max(1L, as.integer(ceiling(min_block / block_size)))

  if (pool == 1L) {
    return(list(
      fold_indices = fold_indices,
      fold_list = sort(unique(fold_indices))
    ))
  }

  if (n_folds %% pool != 0) {
    warning(sprintf(
      paste0("holdout_blocks(): %d folds don't divide into pools of %d - the ",
             "last block will be short. Check V and min_block."),
      n_folds, pool
    ))
  }

  pooled <- as.integer(ceiling(fold_indices / pool))
  message(sprintf(
    "holdout_blocks(): pooling %d folds of ~%.0f into %d blocks of ~%.0f",
    n_folds, block_size, length(unique(pooled)), block_size * pool
  ))

  list(fold_indices = pooled, fold_list = sort(unique(pooled)))
}

#' One stratified train/test split, for the single-split arm
#'
#' Stratified on Y via caret, the same package and the same reasoning as
#' split_folds() - the candidate fits are sensitive enough to the outcome
#' distribution that an unstratified draw would add split-to-split variance
#' this study is not trying to measure.
#'
#' @param Y outcome vector, used only for stratification
#' @param frac fraction of rows to train on, default SPLIT_TRAIN_FRAC
#' @return list of `train` and `test` row indices
single_split <- function(Y, frac = SPLIT_TRAIN_FRAC) {
  train <- caret::createDataPartition(y = Y, p = frac, list = FALSE)[, 1]
  list(train = sort(train), test = sort(setdiff(seq_along(Y), train)))
}
