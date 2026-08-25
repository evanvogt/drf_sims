##########
# title: nuisance-arm pass over existing results - model evaluation study
##########
# A SECOND PASS, not a rerun. me_analysis.R has already produced 358 of 360
# replicates; this reads each one, adds the two new nuisance arms
# (`cv_shared`, `holdout`), carries `whole` and the retired-but-kept
# `cv_indep` through untouched, and writes the result to a parallel tree.
#
# WHY THIS DOES NOT NEED A RERUN. Everything the new arms consume is already
# saved in each res_sim_<run>.RDS (me_analysis.R's results list): data$Y/W/X,
# truth, fold_info, and each candidate's tau. So:
#
#   - the DGM is never re-run. Y/W/X are read, not regenerated. That also
#     sidesteps replaying the RNG stream, which matters because split_folds()
#     consumes it and R/dgm_scenarios.R's header makes draw order a contract.
#   - the 9 candidates are never refit. tau passes through unchanged, so all
#     four arms score the IDENTICAL candidate fits. That is the whole point:
#     the comparison is controlled, differing in the nuisance and nothing
#     else. (me_split.R is the one arm that does refit, by design.)
#   - fold_info is read from disk, never re-derived. Re-deriving it would mean
#     reproducing caret::createFolds() at exactly the right point in the RNG
#     stream, which is a contract this script has no reason to depend on.
#
# WHY A PARALLEL TREE rather than rewriting res_sim_*.RDS in place: those 358
# files are the only copy of a finished run. The per-run objects are plain
# vectors and data.frames (me_collect.R's header), so duplicating them is
# cheap, and it means me_strategies_verify.R can diff old against new rather
# than having to trust that an in-place edit was inert.

library(dplyr)
library(xgboost)
library(h2o)
library(caret)
library(here)

path <- here()

source(here("R", "utils.R"))
source(here("R", "pipeline.R"))
source(here("model_evaluation", "me_utils.R"))
source(here("model_evaluation", "me_nuisance.R"))
source(here("model_evaluation", "me_config.R"))

args <- as.numeric(commandArgs(trailingOnly = TRUE))
i <- args[1]
# no candidate fitting here, so no `workers` knob - n_cores (XGBoost nthread /
# H2O nthreads) is the only one that does anything.
n_cores <- if (length(args) >= 2 && !is.na(args[2])) args[2] else 5
h2o_mem <- "10G"

param <- study$grid[i, ]
print(param)

scenario <- param$scenario
n <- param$n
run <- param$run

# ---- read the source replicate ---------------------------------------------
src_file <- file.path(
  combo_dir(study, param), paste0("res_sim_", run, ".RDS")
)

# The 2 permanently-failed runs have no source file. Exit 0 rather than
# erroring: they are a known, deliberate exclusion, and a non-zero exit would
# put them in this pass's failed_ids list to be resubmitted forever.
if (!file.exists(src_file)) {
  message(
    "no source results at ", src_file,
    " - this is one of the runs excluded from the main study. Skipping."
  )
  quit(save = "no", status = 0)
}

old <- readRDS(src_file)

Y <- old$data$Y
W <- old$data$W
X <- old$data$X

stopifnot(
  "source file has no fold_info" = !is.null(old$fold_info),
  "source file has no nuisances" = !is.null(old$nuisances)
)

# ---- the arms --------------------------------------------------------------
# cv_shared uses the CANDIDATES' fold assignment - the point of the arm is
# that the nuisance for row i trains on exactly the rows the candidate's
# fold-model for row i trained on. holdout inverts that: it sees only the
# fold the candidate held out. Both are derived from fold_info; neither
# draws folds, so neither depends on the RNG.
cand_folds <- old$fold_info
hold_folds <- holdout_blocks(cand_folds$fold_indices)

arms <- list(
  cv_shared = nuisance_arm_spec("cv", cand_folds),
  holdout   = nuisance_arm_spec("holdout", hold_folds)
)

# model_seed is keyed off the ARRAY INDEX, which is unique per
# (scenario, n, run) because the grid rows are - so this pass is reproducible
# without having to reconstruct where in me_analysis.R's stream its own
# model_seed was drawn. It cannot affect the passed-through arms, which are
# copied rather than recomputed.
setup_rng_stream(i)
model_seed <- sample.int(2^31 - 1, 1)

new_nuis <- run_nuisance_arms(
  X, Y, W,
  arms = arms,
  n_cores = n_cores,
  mem = h2o_mem,
  model_seed = model_seed
)

# ---- merge old and new arms ------------------------------------------------
# Ordered by NUISANCE_ARMS so the score columns me_metrics.R emits come out in
# a stable order regardless of which arms were computed where. `cv` is renamed
# to `cv_indep` here and nowhere else - me_analysis.R still writes `cv`, so
# that rerunning one of the two failures produces a file this script can read.
merge_arms <- function(old_pipeline, new_pipeline) {
  merged <- list(
    whole     = old_pipeline$whole,
    cv_indep  = old_pipeline$cv,
    cv_shared = new_pipeline$cv_shared,
    holdout   = new_pipeline$holdout
  )
  stopifnot(
    "merged arm names disagree with me_config.R's NUISANCE_ARMS" =
      identical(names(merged), NUISANCE_ARMS)
  )
  merged
}

nuisances <- lapply(names(old$nuisances), function(pipeline) {
  merge_arms(old$nuisances[[pipeline]], new_nuis[[pipeline]])
})
names(nuisances) <- names(old$nuisances)

# ---- propensity diagnostics ------------------------------------------------
# calculate_pseudos() divides by pi * (1 - pi) with no trimming, and the
# holdout arm fits pi on 25-100 rows, so its predictions sit much closer to
# 0/1 than the whole or cv arms' do - which can make phi's AIPW correction
# blow up and dominate the DR risk. The formula is deliberately NOT changed
# (trimming one arm and not the others would make them non-comparable); this
# records the exposure per arm so the decision can be made from measured
# numbers rather than from the worry.
pi_diagnostics <- function(nuis) {
  bind_rows(lapply(names(nuis), function(pipeline) {
    bind_rows(lapply(names(nuis[[pipeline]]), function(arm) {
      pi <- nuis[[pipeline]][[arm]]$pi
      wt <- 1 / (pi * (1 - pi))
      tibble::tibble(
        pipeline = pipeline, arm = arm,
        pi_min = min(pi), pi_max = max(pi),
        pi_q01 = unname(quantile(pi, 0.01)),
        pi_q99 = unname(quantile(pi, 0.99)),
        weight_max = max(wt), weight_mean = mean(wt),
        n_extreme = sum(pi < 0.05 | pi > 0.95)
      )
    }))
  }))
}

pi_diag <- pi_diagnostics(nuisances)
print(as.data.frame(pi_diag))

# ---- write -----------------------------------------------------------------
# The 9 candidates, data, truth and both fold_infos come straight from `old`,
# so the output is a complete replicate that me_collect.R/me_metrics.R can
# read with no knowledge that it was assembled in two passes.
results <- old
results$nuisances <- nuisances
results$holdout_block_info <- hold_folds
results$pi_diagnostics <- pi_diag
results$strategies_model_seed <- model_seed

output_dir <- combo_dir(study_strat, param)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(output_dir, paste0("res_sim_", run, ".RDS")))

print("Nuisance-arm pass completed!")
