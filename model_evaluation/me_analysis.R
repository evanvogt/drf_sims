##########
# title: script for running the model evaluation study - one replicate per run
##########
# Replaces the old sim_eval.R. Structure mirrors sample_size/correlated/continuous/cts_corr_analysis.R
# (grid comes from the study config, never re-typed here; setup_rng_stream();
# save raw ingredients only, no metrics computed inline - see me_metrics.R).
#
# CLI-argument handling mirrors crossfitting/cf_analysis.R instead: optional
# trailing args for the resource-tuning knobs (workers, n_cores), defaulting
# to sane values so a bare `Rscript me_analysis.R 1` still works as a local
# smoke test. jobscripts/me_1.sh supplies both explicitly.

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

# functions
source(here("R", "utils.R"))
source(here("model_evaluation", "me_dgms.R"))
source(here("model_evaluation", "me_utils.R"))
source(here("model_evaluation", "me_models.R"))
source(here("model_evaluation", "me_nuisance.R"))
source(here("model_evaluation", "me_config.R"))

# simulation parameters
args <- as.numeric(commandArgs(trailingOnly = TRUE))
i <- args[1]
# workers/n_cores default to 4/5 so a bare `Rscript me_analysis.R 1` still
# works as a local smoke test; me_1.sh supplies both explicitly - see the note
# in README.md on sizing the array job.
workers <- if (length(args) >= 2 && !is.na(args[2])) args[2] else 4
n_cores <- if (length(args) >= 3 && !is.na(args[3])) args[3] else 5
h2o_mem <- "10G"

param <- study$grid[i, ]
print(param)

scenario <- param$scenario
n <- param$n
run <- param$run

# candidate models are now single-crossfit (see me_models.R's header), which
# trains each fold's stage-2 model on (V-1)/V of the data rather than
# double-crossfitting's (V-2)/V - so, unlike continuous/, there's no need to
# shrink V at n=250 to keep the training set from getting too small.
n_folds <- 10L

# set up simulation seed
setup_rng_stream(run)

# data generation
gen <- generate_me_scenario_data(scenario, n)
data <- gen$dataset

design <- prepare_design_matrix(data)
Y <- design$Y
W <- design$W
X <- design$X

# k-folds for the candidates' single crossfit. The only fold draw: the scoring
# arms that need folds (cv_shared, holdout) derive them from this one in
# me_strategies.R. A second, independent draw for the scoring nuisance (the
# cv_indep arm) was removed before the 2026-10-08 rerun - see me_config.R's
# NUISANCE_ARMS.
kfolds <- split_folds(Y, k = n_folds)

fold_indices <- kfolds$fold_indices
fold_list <- kfolds$fold_list

# candidate models
metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)

model_list <- run_all_candidate_models(Y, W, X, fold_indices, fold_list)

# nuisance-evaluation pipelines, used only to score the candidates above.
# Only the `whole` arm here; me_strategies.R adds cv_shared and holdout.
# model_seed is drawn AFTER data generation, the fold draw, and
# candidate-model fitting - keep it in this exact position, since draw order
# from setup_rng_stream() is part of the reproducibility contract (see
# R/dgm_scenarios.R's header for the same principle applied to the DGM
# itself).
model_seed <- sample.int(2^31 - 1, 1)

nuisances <- run_nuisance_arms(
  X, Y, W,
  arms = list(whole = nuisance_arm_spec("whole")),
  n_cores = n_cores, mem = h2o_mem, model_seed = model_seed
)

# save fitted models, nuisances, data and truth. The 9 candidate models are
# spliced into the TOP LEVEL of results (not nested under a "models" key) -
# me_metrics.R's use of R/metrics.R::compute_metrics() depends on this, see
# its header comment.
results <- c(model_list, list(
  data = list(Y = Y, W = W, X = X),
  truth = gen$truth,
  fold_info = kfolds,
  nuisances = nuisances
))

# Save results - to combo_dir(), the path me_strategies.R and me_split.R read
# their source runs from, so the three trees follow study$res_path together
output_dir <- combo_dir(study, param)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(output_dir, paste0("res_sim_", run, ".RDS")))

print("Simulation completed!")
