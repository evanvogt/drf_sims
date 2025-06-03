collate_predictions <- function(fold_list, fold_pairs, fold_indices, cross_fits, target) {
  lapply(fold_list, function(fold) {
    predictions <- rep(NA, length.out = length(fold_indices))
    for (j in seq_along(fold_pairs)) {
      fold_pair <- fold_pairs[[j]]
      if (fold %in% fold_pair) next
      
      predictions[fold_indices %in% fold_pair] <- cross_fits[[j]][[target]]
    }
    predictions[fold_indices == fold] <- NA
    predictions
  }) %>% simplify2array()
}