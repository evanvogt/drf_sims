##########
# title: shared low-level helpers
##########
# Sourced by every study, directly. The repo-root utils.R is a shim onto this
# file, kept only because case_study/ relies on whatever drives it to have
# loaded that shim first - case_study/ is outside the scope of the restructure.

require(parallel)

# ---- reproducibility --------------------------------------------------------

#' Give each simulation run an independent RNG stream
#'
#' Uses L'Ecuyer-CMRG substreams so that run `i` is reproducible on its own,
#' independently of how many workers the job happens to use.
#'
#' @param trial_num run index; advances the stream this many times
#' @param seed base seed shared by every run
setup_rng_stream <- function(trial_num, seed = 1998) {
  RNGkind("L'Ecuyer-CMRG")

  set.seed(seed)

  for (i in 1:trial_num) {
    .Random.seed <<- parallel::nextRNGStream(.Random.seed)
  }
}

# ---- fold bookkeeping -------------------------------------------------------

#' Assemble double-crossfitting predictions into an n x V matrix
#'
#' Column `k` holds predictions from the nuisance models that never saw fold `k`,
#' and is NA on fold `k`'s own rows. The stage-2 model for fold `k` therefore
#' trains on `po_matrix[, k]` without any leakage from the fold it will predict.
#'
#' @param fold_list unique fold labels
#' @param fold_pairs list of fold pairs, as from utils::combn(fold_list, 2)
#' @param fold_indices length-n vector of fold membership
#' @param reslist per-fold-pair results, in the same order as fold_pairs
#' @param target name of the field to pull out of each element of reslist
collate_predictions <- function(fold_list, fold_pairs, fold_indices, reslist, target) {
  simplify2array(lapply(fold_list, function(fold) {
    predictions <- rep(NA, length(fold_indices))
    for (j in seq_along(fold_pairs)) {
      if (fold %in% fold_pairs[[j]]) {
        predictions[fold_indices %in% fold_pairs[[j]]] <- reslist[[j]][[target]]
      }
    }
    predictions[fold_indices == fold] <- NA
    predictions
  }))
}

#' Reassemble leave-one-fold-out results into full-length vectors
#'
#' The single-crossfitting counterpart of collate_predictions: each element of
#' reslist carries a `fold` label and covers exactly that fold's rows.
#'
#' @param reslist per-fold results, each with a `fold` element
#' @param fold_indices length-n vector of fold membership
#' @param targets names of the fields to extract
#' @return named list of length-n numeric vectors
scatter_folds <- function(reslist, fold_indices, targets) {
  out <- lapply(targets, function(nm) {
    v <- rep(NA_real_, length(fold_indices))
    for (res in reslist) v[fold_indices == res$fold] <- res[[nm]]
    v
  })
  names(out) <- targets
  out
}

# ---- misc -------------------------------------------------------------------

# clamp propensities away from 0 and 1 so the DR denominator cannot explode
trim_ps <- function(p, lo = 0.05, hi = 0.95) pmin(pmax(p, lo), hi)

# evaluate an expression, returning its value alongside elapsed seconds
timed <- function(expr) {
  t <- system.time(val <- expr)
  list(value = val, time = unname(t["elapsed"]))
}
