##########
# Title: diagnose the SuperLearner "missing data" failure in the split
# pseudo-observation DR-learner (nuisance_pseudo_sl_split / stage_2_sl_vec)
##########
# Context: all_cate_surv_models() aborts with
#   SuperLearner(): missing data is currently not supported
# when it reaches results$sl_dr_split (surv_models.R:307-339). The standard
# cross-fit pseudo-value paths (pseudo_crossfit, pseudo_double_crossfit)
# guard against pseudoyl()/pseudomean() returning NA by falling back to
# pseudo_whole; the split-pseudo-obs path (nuisance_pseudo_sl_split,
# compute_split_pseudoyl, compute_split_pseudomean) has no such guard
# anywhere. This script isolates exactly which unguarded vector first
# carries NAs, at what rate, and why - it does not fix anything.
#
# ROOT CAUSE (confirmed by stepping into pseudo:::ci.omit, the internal
# leave-one-out competing-risks engine behind pseudoyl()'s cause-specific
# pseudo-values - pseudomean()'s single-cause surv.omit() does not have
# this bug, which is why RMSTc never shows NAs below):
#
#   ci.omit() computes, for every individual i in one shot via a risk-set
#   matrix trick, the leave-i-out cumulative-incidence curve. For row
#   i = the individual with the LARGEST observed time in the sample, the
#   leave-one-out risk set can hit exactly 0 at the tail (nobody else is
#   still at risk once that individual's own time is excluded), producing
#   a 0/0 = NaN both in the shared KM-survival term `km` and in the
#   cause-specific numerator ratio `cumi <- N1/Y`. The package patches the
#   *km* matrix's trailing NaNs by forward-filling (holding the survival
#   curve flat past the last supported time) - but never patches the
#   matching NaNs in `cumi`. Since the cause-specific CI is
#   `cumsum(cumi * km)` along that row, the single unpatched NaN poisons
#   every subsequent cumsum entry, and pseudoyl() returns NA for exactly
#   that one individual: whoever holds the maximum observed time in the
#   sample being decomposed.
#
#   compute_split_pseudoyl(Y_km, D_km, Y_val, D_val, horizon) evaluates
#   pseudoyl(c(Y_km, Y_val[j]), ...) and keeps only the *last* (the
#   appended validation individual's) pseudo-value. So every time
#   Y_val[j] exceeds max(Y_km), the appended individual becomes exactly
#   this "row howmany" case and its pseudo-value comes back NA - this is
#   exactly the pattern in step 2/3 below (offending Y always > the KM
#   set's max Y). The same failure mode already affects the *training-fold*
#   pseudoyl() calls too (see the nonzero na_train_* columns for
#   seed=2 in step 4) for the same reason - whoever has the max Y within
#   that training fold - which is precisely why pseudo_crossfit/
#   pseudo_double_crossfit already carry an is.na()->pseudo_whole
#   fallback (surv_models.R:601-633, 692-713): this is a known, structural
#   limitation of the `pseudo` package's ci.omit(), not a data/sampling bug.
#
# Run from competing_risk/:
#   Rscript surv_dr_split_na_diagnose.R

library(dplyr)
library(furrr)
library(pseudo)
library(SuperLearner)

source("surv_dgm.R")
source("surv_models.R")

cat(sprintf(
  "pseudo package version: %s\nSuperLearner package version: %s\nR version: %s\n\n",
  as.character(packageVersion("pseudo")),
  as.character(packageVersion("SuperLearner")),
  R.version.string
))

# production settings (competing_risk/surv_analysis.R, README.md)
horizon <- 28
n <- 500
n_folds <- 10

plan(sequential) # keep this diagnostic single-threaded and deterministic

# =============================================================================
cat("=== 1. Reproduce the failure once ===\n")

set.seed(1)
gen <- generate_surv_data(scenario = 6, n = n, censoring = FALSE)
data <- gen$dataset

X <- as.matrix(data[, !(names(data) %in% c("Y", "D", "W"))])
Y <- data$Y
D <- data$D
W <- data$W

fold_indices <- sort(seq(nrow(X)) %% n_folds) + 1
fold_list <- unique(fold_indices)

nuisances <- tryCatch(
  {
    nuisance_pseudo_sl_split(
      X,
      Y,
      D,
      W,
      horizon,
      fold_indices,
      fold_list,
      DEFAULT_SL_LIBRARY
    )
  },
  error = function(e) {
    cat(sprintf("  nuisance_pseudo_sl_split() errored: %s\n\n", conditionMessage(e)))
    NULL
  }
)

if (!is.null(nuisances)) {
  cat(sprintf(
    "  nuisance_pseudo_sl_split() completed. NAs in po_RMTL1/po_RMTL2/po_RMSTc: %d / %d / %d\n",
    sum(is.na(nuisances$po_RMTL1)),
    sum(is.na(nuisances$po_RMTL2)),
    sum(is.na(nuisances$po_RMSTc))
  ))
  cat(
    "  (the NAs seen in step 2's split_RMTL* enter po via the unguarded formula\n",
    "  at surv_models.R:1331-1333, not as a direct SuperLearner() Y argument -\n",
    "  so nuisance_pseudo_sl_split() itself does not error even when it should.)\n\n"
  )

  # This is the actual downstream call site the user reported (surv_models.R:319-339).
  stage2 <- tryCatch(
    {
      stage_2_sl_vec(
        X,
        nuisances$po_RMTL1,
        fold_indices,
        fold_list,
        DEFAULT_SL_LIBRARY
      )
      "no error"
    },
    error = function(e) conditionMessage(e)
  )
  cat(sprintf("  stage_2_sl_vec(po_RMTL1) result: %s\n\n", stage2))
}

# =============================================================================
cat("=== 2. Per-fold isolation of the four unguarded quantities ===\n")
# Replicates the exact split logic from nuisance_pseudo_sl_split
# (surv_models.R:1272-1305) but stops short of calling SuperLearner, so we
# can is.na()-count each intermediate quantity independently.

D_int <- as.integer(D)
Dc <- as.integer(D %in% c(1, 2))

na_table <- data.frame()
offenders <- list()

for (i in seq_along(fold_list)) {
  fold <- fold_list[i]
  km_fold <- fold_list[(i %% n_folds) + 1]

  in_val <- fold_indices == fold
  in_km <- fold_indices == km_fold
  in_train <- !in_val & !in_km

  ps_RMTL_train <- pseudoyl(Y[in_train], D_int[in_train], horizon)
  ps_RMSTc_train <- pseudomean(Y[in_train], Dc[in_train], horizon)

  split_RMTL <- compute_split_pseudoyl(
    Y[in_km],
    D_int[in_km],
    Y[in_val],
    D_int[in_val],
    horizon
  )
  split_RMSTc <- compute_split_pseudomean(
    Y[in_km],
    Dc[in_km],
    Y[in_val],
    Dc[in_val],
    horizon
  )

  row <- data.frame(
    fold = fold,
    n_train = sum(in_train),
    n_km = sum(in_km),
    n_val = sum(in_val),
    na_train_RMTL1 = sum(is.na(ps_RMTL_train$pseudo$cause1)),
    na_train_RMTL2 = sum(is.na(ps_RMTL_train$pseudo$cause2)),
    na_train_RMSTc = sum(is.na(ps_RMSTc_train)),
    na_split_RMTL1 = sum(is.na(split_RMTL$cause1)),
    na_split_RMTL2 = sum(is.na(split_RMTL$cause2)),
    na_split_RMSTc = sum(is.na(split_RMSTc))
  )
  na_table <- rbind(na_table, row)

  # stash offending rows for step 3
  val_idx <- which(in_val)
  train_idx <- which(in_train)
  km_idx <- which(in_km)

  if (any(is.na(ps_RMTL_train$pseudo$cause1))) {
    bad <- train_idx[is.na(ps_RMTL_train$pseudo$cause1)]
    offenders[[length(offenders) + 1]] <- list(
      fold = fold,
      quantity = "train_RMTL1",
      idx = bad,
      Y = Y[bad],
      D = D[bad]
    )
  }
  if (any(is.na(ps_RMTL_train$pseudo$cause2))) {
    bad <- train_idx[is.na(ps_RMTL_train$pseudo$cause2)]
    offenders[[length(offenders) + 1]] <- list(
      fold = fold,
      quantity = "train_RMTL2",
      idx = bad,
      Y = Y[bad],
      D = D[bad]
    )
  }
  if (any(is.na(ps_RMSTc_train))) {
    bad <- train_idx[is.na(ps_RMSTc_train)]
    offenders[[length(offenders) + 1]] <- list(
      fold = fold,
      quantity = "train_RMSTc",
      idx = bad,
      Y = Y[bad],
      D = D[bad]
    )
  }
  if (any(is.na(split_RMTL$cause1))) {
    bad <- val_idx[is.na(split_RMTL$cause1)]
    offenders[[length(offenders) + 1]] <- list(
      fold = fold,
      quantity = "split_RMTL1",
      idx = bad,
      Y = Y[bad],
      D = D[bad],
      km_Y_range = range(Y[km_idx]),
      km_events = table(D[km_idx])
    )
  }
  if (any(is.na(split_RMTL$cause2))) {
    bad <- val_idx[is.na(split_RMTL$cause2)]
    offenders[[length(offenders) + 1]] <- list(
      fold = fold,
      quantity = "split_RMTL2",
      idx = bad,
      Y = Y[bad],
      D = D[bad],
      km_Y_range = range(Y[km_idx]),
      km_events = table(D[km_idx])
    )
  }
  if (any(is.na(split_RMSTc))) {
    bad <- val_idx[is.na(split_RMSTc)]
    offenders[[length(offenders) + 1]] <- list(
      fold = fold,
      quantity = "split_RMSTc",
      idx = bad,
      Y = Y[bad],
      D = D[bad],
      km_Y_range = range(Y[km_idx]),
      km_events = table(D[km_idx])
    )
  }
}

print(na_table, row.names = FALSE)
cat("\n")

# =============================================================================
cat("=== 3. Characterize the NA-producing observations ===\n")

if (length(offenders) == 0) {
  cat(
    "  No NAs found in scenario 6 / n=500 / censoring=FALSE / seed=1.\n",
    "  (Failure in step 1, if any, must come from elsewhere - re-check.)\n\n"
  )
} else {
  for (o in offenders) {
    cat(sprintf(
      "  fold %d, %s: %d NA(s) at index/indices %s\n",
      o$fold,
      o$quantity,
      length(o$idx),
      paste(o$idx, collapse = ", ")
    ))
    cat(sprintf(
      "    offending Y: %s\n    offending D: %s\n",
      paste(round(o$Y, 2), collapse = ", "),
      paste(o$D, collapse = ", ")
    ))
    if (!is.null(o$km_Y_range)) {
      cat(sprintf(
        "    KM-set Y range: [%.2f, %.2f], horizon = %d\n",
        o$km_Y_range[1],
        o$km_Y_range[2],
        horizon
      ))
      cat("    KM-set D counts: ")
      print(o$km_events)
    }
    cat("\n")
  }
}

# =============================================================================
cat("=== 4. Breadth check: scenario x censoring x seed ===\n")

breadth <- data.frame()
scenarios <- c(1, 6)
censoring_opts <- c(FALSE, TRUE)
seeds <- c(1, 2)

for (sc in scenarios) {
  for (cens in censoring_opts) {
    for (sd in seeds) {
      set.seed(sd)
      g <- generate_surv_data(scenario = sc, n = n, censoring = cens)
      d <- g$dataset
      Yb <- d$Y
      Db <- d$D
      D_int_b <- as.integer(Db)
      Dc_b <- as.integer(Db %in% c(1, 2))

      fi <- sort(seq(nrow(d)) %% n_folds) + 1
      fl <- unique(fi)

      total_na_train <- 0
      total_na_split <- 0

      for (i in seq_along(fl)) {
        fold <- fl[i]
        km_fold <- fl[(i %% n_folds) + 1]
        in_val <- fi == fold
        in_km <- fi == km_fold
        in_train <- !in_val & !in_km

        ps_t <- pseudoyl(Yb[in_train], D_int_b[in_train], horizon)
        ps_c <- pseudomean(Yb[in_train], Dc_b[in_train], horizon)
        total_na_train <- total_na_train +
          sum(is.na(ps_t$pseudo$cause1)) +
          sum(is.na(ps_t$pseudo$cause2)) +
          sum(is.na(ps_c))

        sp_t <- compute_split_pseudoyl(
          Yb[in_km],
          D_int_b[in_km],
          Yb[in_val],
          D_int_b[in_val],
          horizon
        )
        sp_c <- compute_split_pseudomean(
          Yb[in_km],
          Dc_b[in_km],
          Yb[in_val],
          Dc_b[in_val],
          horizon
        )
        total_na_split <- total_na_split +
          sum(is.na(sp_t$cause1)) +
          sum(is.na(sp_t$cause2)) +
          sum(is.na(sp_c))
      }

      breadth <- rbind(
        breadth,
        data.frame(
          scenario = sc,
          censoring = cens,
          seed = sd,
          na_train_total = total_na_train,
          na_split_total = total_na_split
        )
      )
    }
  }
}

print(breadth, row.names = FALSE)
cat("\n")

# =============================================================================
cat("=== 5. Summary ===\n")

any_train_na <- sum(na_table$na_train_RMTL1) +
  sum(na_table$na_train_RMTL2) +
  sum(na_table$na_train_RMSTc)
any_split_na <- sum(na_table$na_split_RMTL1) +
  sum(na_table$na_split_RMTL2) +
  sum(na_table$na_split_RMSTc)

cat(sprintf(
  paste0(
    "  scenario 6 / n=500 / 10 folds / horizon=28 / censoring=FALSE / seed=1:\n",
    "    training-fold pseudoyl/pseudomean (surv_models.R:1288-1289) NAs: %d\n",
    "    split pseudoyl/pseudomean (surv_models.R:1292-1305, via compute_split_*",
    " at 929-965) NAs: %d\n\n",
    "  breadth check (scenarios 1/6, censoring T/F, 2 seeds) total NAs:\n",
    "    training-fold: %d   split: %d\n\n",
    "  -> whichever total is nonzero pinpoints the call site: if training-fold\n",
    "     is nonzero, surv_models.R:1288-1289 needs the same pseudo_whole\n",
    "     fallback pattern already used in pseudo_crossfit (lines 601-633); if\n",
    "     split is nonzero, compute_split_pseudoyl/compute_split_pseudomean\n",
    "     (lines 929-965) and their caller (lines 1292-1305) need an analogous\n",
    "     guard with no pseudo_whole equivalent currently available for the\n",
    "     split construction, so a design decision is needed for what to fall\n",
    "     back to (candidates: pseudo_whole itself, or drop the observation\n",
    "     from that fold's stage-2 fit).\n"
  ),
  any_train_na,
  any_split_na,
  sum(breadth$na_train_total),
  sum(breadth$na_split_total)
))
