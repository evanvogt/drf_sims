##########
# title: script for running the model evaluation study - one replicate per run
##########
# Replaces the old sim_eval.R. Structure mirrors continuous/cts_analysis.R
# (grid comes from the study config, never re-typed here; setup_rng_stream();
# save raw ingredients only, no metrics computed inline - see me_metrics.R).
#
# CLI-argument handling mirrors crossfitting/cf_analysis.R instead: optional
# trailing args for the resource-tuning knobs (workers, n_cores), defaulting
# to sane values so a bare `Rscript me_analysis.R 1` still works as a local
# smoke test. jobscripts/me_1.sh supplies real, measured values once
# me_profile_summary.R has run - see me_profile.R/me_profile_summary.R.

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

# path
path <- here()

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
# workers/n_cores default to the pre-profiling values so a bare
# `Rscript me_analysis.R 1` still works as a local smoke test; me_1.sh always
# supplies both once me_profile_summary.R has written them in - see
# me_profile_summary.R and the note in README.md on sizing the array job.
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

# k-folds for cross-fitting/validation. Two independent draws, both from
# split_folds() and both consuming the RNG stream (see model_seed's comment
# below for why draw order matters here): one for the 9 candidate models,
# one for me_nuisance.R's scoring pipelines. Under the old double-crossfit
# candidates the two could safely share one fold_indices - the scoring
# pipeline's leave-one-fold-out nuisance for row i and a double-crossfit
# candidate's tau_hat at row i were never fit on the identical training set.
# Single crossfitting removes that separation: sharing folds would fit
# tau_hat_i and the scoring nuisance at row i on the identical row set,
# correlating their errors and biasing the proxy score toward whatever
# tau_hat the shared split happens to favour. An independent draw keeps the
# scoring honest.
kfolds <- split_folds(Y, k = n_folds)
nuis_folds <- split_folds(Y, k = n_folds)

fold_indices <- kfolds$fold_indices
fold_list <- kfolds$fold_list

# candidate models
metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)

model_list <- run_all_candidate_models(Y, W, X, fold_indices, fold_list)

# nuisance-evaluation pipelines, used only to score the candidates above.
# model_seed is drawn AFTER data generation, both fold draws, and
# candidate-model fitting - keep it in this exact position, since draw order
# from setup_rng_stream() is part of the reproducibility contract (see
# R/dgm_scenarios.R's header for the same principle applied to the DGM
# itself).
model_seed <- sample.int(2^31 - 1, 1)

nuisances <- run_all_nuisance_pipelines(
  X, Y, W, nuis_folds$fold_indices, nuis_folds$fold_list,
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
  nuisance_fold_info = nuis_folds,
  nuisances = nuisances
))

# Save results
output_dir <- file.path(dirname(path), "results", "model_evaluation",
                        paste0("scenario_", scenario), n)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(output_dir, paste0("res_sim_", run, ".RDS")))

print("Simulation completed!")
