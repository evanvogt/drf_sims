###############
# Simulation steps all in one script
###############

# libraries -----
library(here)
library(benchtm)
library(tidyverse)
library(caret)

# paths -----
sim_path <- here("src", "selection_simulation")
utils_path <- here("src", "utils")

# functions ----
source(here(utils_path, "sim_utils.R"))
source(here(utils_path, "model_utils.R"))
source(here(utils_path, "nuisance_utils.R"))
source(here(utils_path, "metric_utils.R"))


# simulation parameters
args <- commandArgs(trailingOnly = TRUE) # if not running on HPC, specify parameters manually
scen_index <- as.numeric(args[1])
sim_index <- as.numeric(args[2])

k <- 10 # number of folds for crossfitting
n_cores <- 5 # option to parallelise with more cores (defaults in functions are single core)
mem <- "10G" # memory allocation for h2o cluster


# parallelisation option
oplan <- plan(multisession, workers = 2)
on.exit(oplan, add = TRUE)
# what are the steps?

# main -----
# 0. simulation seed -----
setup_rng_seed(sim_index)

# 1. data generation -----
data(scen_param)
gen <- generate_data(scen_param, scen_index)
adat <- gen$adat
true_tau <- gen$true_tau

design <- prepare_design_matrix(adat)
Y <- design$Y
W <- design$W
X <- design$X

nvar <- ncol(X)
N <- nrow(X)

# k-folds for cross-fitting/validation
kfolds <- split_folds(Y, k = k)

fold_indices <- kfolds$fold_indices
fold_list <- kfolds$fold_list
fold_pairs <- kfolds$fold_pairs


# 2. candidate models -----
simple <- list(
  hy_rf1 = list(), # Defaults
  hy_rf2 = create_rf_hyperparams(mtry = 30, max.depth = 5),
  hy_rf3 = create_rf_hyperparams(mtry = 10, max.depth = 3),

  hy_net1 = create_net_hyperparams(alpha = 1), # Default Lasso
  hy_net2 = create_net_hyperparams(alpha = 0), # Ridge
  hy_net3 = create_net_hyperparams(alpha = 0.5), # Elastic Net

  hy_SL1 = create_SL_hyperparams(
    method = "method.CC_LS",
    SL.library = list("SL.glmnet", "SL.ranger", "SL.earth", "SL.gam", "SL.mean")
  ),
  hy_SL2 = create_SL_hyperparams(
    method = "method.CC_LS",
    SL.library = list(
      "SL.glmnet",
      "SL.xgboost",
      "SL.cforest",
      "SL.earth",
      "SL.gam",
      "SL.mean"
    )
  ),
  hy_SL3 = create_SL_hyperparams(
    method = "method.CC_LS",
    SL.library = list("SL.svm", "SL.nnet", "SL.mean")
  )
)
list2env(simple, envir = .GlobalEnv)
# 3. model fitting -----
rf_params <- grep("hy_rf", names(.GlobalEnv), value = TRUE)
for (params in rf_params) {
  res_name <- sub("hy_", "fit_", params)
  fitted_model <- fit_rf(
    Y,
    W,
    X,
    hyper_list = get(params),
    fold_indices,
    fold_list,
    fold_pairs
  )
  assign(res_name, fitted_model, envir = .GlobalEnv)
}

net_params <- grep("hy_net", names(.GlobalEnv), value = TRUE)
for (params in net_params) {
  res_name <- sub("hy_", "fit_", params)
  fitted_model <- fit_glmnet(
    Y,
    W,
    X,
    hyper_list = get(params),
    fold_indices,
    fold_list,
    fold_pairs
  )
  assign(res_name, fitted_model, envir = .GlobalEnv)
}

SL_params <- grep("hy_SL", names(simple), value = TRUE)
for (params in SL_params) {
  res_name <- sub("hy_", "fit_", params)
  fitted_model <- fit_SL(
    Y,
    W,
    X,
    hyper_list = get(params),
    fold_indices,
    fold_list,
    fold_pairs
  )
  assign(res_name, fitted_model, envir = .GlobalEnv)
}

# save fitted models, data, truth and seed
models <- grep("^fit_.*\\d$", names(.GlobalEnv), value = TRUE)
model_list <- mget(models, envir = .GlobalEnv)

# 4. estimate nuisance functions for evaluation -----
model_seed <- sample.int(2^31 - 1, 1)

nuisance_models <- list(
  xgb = run_all_xgb_nuisance(X, Y, W, fold_indices, fold_list),
  automl = run_all_automl_nuisance(
    X,
    Y,
    W,
    fold_indices,
    fold_list,
    n_cores = n_cores,
    mem = mem,
    model_seed = model_seed
  )
)

# 5. get evaluation metrics (*2) on all models -----
metrics <- calc_metrics(
  models_list,
  nuisance_models,
  Y,
  W,
  truth_avail = TRUE,
  true_tau
)


# 6. save results -----
results <- list(
  data = list(Y = Y, W = W, X = X, truth = true_tau),
  fold_info = kfolds,
  models = model_list,
  nuisances = nuisance_models,
  metrics = metrics
)

out_dir <- here("src", "selection_simulation", "results")

if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

res_name <- paste0("results_scenario_", scen_index, "_sim_", sim_index, ".RDS")
saveRDS(results, here(out_dir, res_name))
