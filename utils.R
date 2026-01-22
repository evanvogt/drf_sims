require(parallel)

# function for making setting independent simulation streams
setup_rng_stream <- function(trial_num, seed = 1998) {
  RNGkind("L'Ecuyer-CMRG")
  
  set.seed(seed)

  for(i in 1:trial_num) {
    .Random.seed <<- parallel::nextRNGStream(.Random.seed)
  }
}


# collecting up predictions for double crossfitting procedure
collate_predictions <- function(fold_list, fold_pairs, fold_indices, reslist, target) {
  lapply(fold_list, function(fold) {
    predictions <- rep(NA, length(fold_indices))
    for (j in seq_along(fold_pairs)) {
      if (fold %in% fold_pairs[[j]]) {
        predictions[fold_indices %in% fold_pairs[[j]]] <- reslist[[j]][[target]]
      }
    }
    predictions[fold_indices == fold] <- NA
    predictions
  }) %>% simplify2array()
}