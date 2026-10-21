##########
# title: the 80:20 single-split arm - model evaluation study
##########
# The 80:20 arm: the 9 candidates are refit on 80% and the nuisance sees only
# the other 20%, so neither borrows from the other's rows. Every other arm
# scores me_analysis.R's crossfit fits, which have seen all n rows.
#
# A separate tree, not a column beside the other arms: tau_hat, the evaluation
# rows and true_pehe all live on the 20%, so one table would compare scores over
# different rows against differently computed truths. n = 250 is excluded: a
# 50-row evaluation set is too thin to rank 9 models.
#
# `data` and `truth` are stored already restricted to the 20%, so
# compute_metrics() and me_per_model() run unchanged; the full vectors and
# split indices are kept under separate names.

library(dplyr)
library(future)
library(future.apply)
library(ranger)
library(glmnet)
library(SuperLearner)
library(xgboost)
library(h2o)
library(caret)
library(here)

path <- here()

source(here("R", "utils.R"))
source(here("R", "pipeline.R"))
source(here("model_evaluation", "me_utils.R"))
source(here("model_evaluation", "me_models.R"))
source(here("model_evaluation", "me_nuisance.R"))
source(here("model_evaluation", "me_config.R"))

args <- as.numeric(commandArgs(trailingOnly = TRUE))
i <- args[1]
workers <- if (length(args) >= 2 && !is.na(args[2])) args[2] else 4
n_cores <- if (length(args) >= 3 && !is.na(args[3])) args[3] else 5
h2o_mem <- "10G"

param <- study_split$grid[i, ]
print(param)

scenario <- param$scenario
n <- param$n
run <- param$run

# folds for the crossfit within the 80%. Still 10: the training set here is
# 0.8n, so at the smallest n this arm runs (500) a fold is 40 rows and the
# stage-1 models train on 360 - no reason to shrink V, for the same reason
# me_analysis.R doesn't.
n_folds <- 10L

# ---- read the source replicate ---------------------------------------------
# Same reasoning as me_strategies.R: the data is read, never regenerated, so
# this arm runs against byte-identical Y/W/X to every other arm. Only the
# split and the fitting differ.
src_file <- file.path(
  combo_dir(study, param), paste0("res_sim_", run, ".RDS")
)

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
true_tau <- old$truth$tau

stopifnot(
  "source file has no truth$tau" = !is.null(true_tau),
  "truth$tau is not length n" = length(true_tau) == length(Y)
)

# ---- the split -------------------------------------------------------------
# Keyed off the array index, which is unique per (scenario, n, run) because
# study_split's grid rows are - so the split is reproducible without
# depending on where in me_analysis.R's stream anything was drawn. Saved
# below regardless, so nothing downstream has to re-derive it.
setup_rng_stream(i)
sp <- single_split(Y)
train <- sp$train
test <- sp$test

cat(sprintf(
  "single split: %d train / %d test (%.0f%% / %.0f%%)\n",
  length(train), length(test),
  100 * length(train) / length(Y), 100 * length(test) / length(Y)
))

# ---- candidates, fit on the 80% only ---------------------------------------
metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)

model_list <- run_all_candidate_models_split(
  Y, W, X, train, test, n_folds = n_folds
)

stopifnot(
  "a candidate returned the wrong number of predictions" =
    all(vapply(model_list, function(m) length(m$tau) == length(test), logical(1)))
)

# ---- the nuisance, fit on the 20% only -------------------------------------
# `whole` applied to the test rows and nothing else. Fit and predicted on the
# same rows, so this arm is not row-level honest - the same trade the
# `holdout` arm makes, and for the same reason: with the candidate excluded
# from these rows entirely, there is no larger set left to fit an honest
# nuisance from. cv_shared in the other tree is the row-honest comparator.
model_seed <- sample.int(2^31 - 1, 1)

nuisances <- run_nuisance_arms(
  X[test, , drop = FALSE], Y[test], W[test],
  arms = list(split = nuisance_arm_spec("whole")),
  n_cores = n_cores,
  mem = h2o_mem,
  model_seed = model_seed
)

# ---- assemble --------------------------------------------------------------
results <- c(model_list, list(
  # restricted to the evaluation rows - see the header note on metrics
  data = list(Y = Y[test], W = W[test], X = X[test, , drop = FALSE]),
  truth = list(tau = true_tau[test]),
  nuisances = nuisances,
  # Traceability, under names compute_metrics() ignores (it takes only
  # intersect(names(sim_res), models)). The row indices are enough to
  # reconstruct everything against the source tree, so the full Y/W/X are
  # not duplicated here - they already exist, unchanged, in
  # results/model_evaluation, and this arm reads rather than regenerates them.
  split_info = list(
    train = train, test = test, n_folds = n_folds,
    source_file = src_file
  ),
  split_model_seed = model_seed
))

output_dir <- combo_dir(study_split, param)
save_run(results, output_dir, run)

print("Single-split arm completed!")
