##########
# title: nuisance-arm pass over existing results - model evaluation study
##########
# A second pass over me_analysis.R's replicates: reads each one, fits the three
# nuisance arms (`whole`, `cv_shared`, `holdout`) and writes a parallel tree.
# Only the main tree's candidates, data, truth and fold_info are read; `whole`
# is recomputed so all arms come from the current pipelines (DEVLOG bug R).
#
# A second pass because nothing is regenerated: no DGM re-run or fold re-draw
# (so no RNG replay) and no candidate refit, so all three arms score the
# identical fits. A parallel tree because res_sim_*.RDS is the only copy of a
# finished run, and me_strategies_verify.R can diff old against new.

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
  "source file has no fold_info" = !is.null(old$fold_info)
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
  whole     = nuisance_arm_spec("whole"),
  cv_shared = nuisance_arm_spec("cv", cand_folds),
  holdout   = nuisance_arm_spec("holdout", hold_folds)
)

# model_seed is keyed off the ARRAY INDEX, which is unique per
# (scenario, n, run) because the grid rows are - so this pass is reproducible
# without having to reconstruct where in me_analysis.R's stream its own
# model_seed was drawn.
setup_rng_stream(i)
model_seed <- sample.int(2^31 - 1, 1)

new_nuis <- run_nuisance_arms(
  X, Y, W,
  arms = arms,
  n_cores = n_cores,
  mem = h2o_mem,
  model_seed = model_seed
)

# ---- order the arms --------------------------------------------------------
# Ordered by NUISANCE_ARMS so the score columns me_metrics.R emits come out in
# a stable order.
nuisances <- lapply(new_nuis, function(pipeline) {
  ordered <- pipeline[NUISANCE_ARMS]
  stopifnot(
    "arm names disagree with me_config.R's NUISANCE_ARMS" =
      identical(names(ordered), NUISANCE_ARMS)
  )
  ordered
})

# ---- propensity diagnostics ------------------------------------------------
# On the RAW pi, before calculate_pseudos()'s trimming to [0.05, 0.95]. The
# holdout arm fits pi on 25-100 rows and predicts the rows it was fit on, so
# its predictions sit much closer to 0/1 than the cv_shared arm's do;
# n_extreme is how many rows per arm the trimming moved. Untrimmed, exact 0/1
# predictions made phi undefined (DEVLOG bug S).
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
# The 9 candidates, data, truth and fold_info come straight from `old`,
# so the output is a complete replicate that me_collect.R/me_metrics.R can
# read with no knowledge that it was assembled in two passes.
results <- old
results$nuisances <- nuisances
results$holdout_block_info <- hold_folds
results$pi_diagnostics <- pi_diag
results$strategies_model_seed <- model_seed

output_dir <- combo_dir(study_strat, param)
save_run(results, output_dir, run)

print("Nuisance-arm pass completed!")
