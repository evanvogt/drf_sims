##########
# title: optimal sample.fraction calibration simulation - continuous outcome
##########

# libraries
library(here)
library(dplyr)

# path
path <- here()

# functions
source(here("utils.R"))
source(here("confidence_intervals", "continuous", "cts_ci_dgms.R"))
source(here("confidence_intervals", "continuous", "cts_ci_models.R"))
source(here("confidence_intervals", "optimal_sf", "cts_ci_sf_calibration.R"))

# simulation parameters
i <- as.numeric(commandArgs(trailingOnly = T))

n_folds <- 10
CI_boot <- 200
alpha   <- 0.05
workers <- 2

# no CI_sf axis — that is what we are finding
params <- expand.grid(
  scenario = c(1:10),
  n        = c(500, 1000),
  run      = c(1:100),
  stringsAsFactors = FALSE
)

param    <- params[i, ]
print(param)

scenario <- param$scenario
n        <- param$n
run      <- param$run

# set up simulation seed
setup_rng_stream(run)

# data generation
gen  <- generate_continuous_scenario_data(scenario, n)
data <- gen$dataset

X            <- as.matrix(data[, -c(1:2)])
Y            <- data$Y
W            <- data$W
n_obs        <- nrow(X)
fold_indices <- sort(seq(n_obs) %% n_folds) + 1
fold_list    <- unique(fold_indices)
fold_pairs   <- utils::combn(fold_list, 2, simplify = FALSE)

# set up parallelisation
metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)

# step 1: lightweight fit — nuisances + tau.hat only, no bootstrap CIs
cat("Computing nuisance functions...\n")
nuisances_rf <- nuisance_rf(X, Y, W, fold_indices, fold_pairs)

cat("Estimating tau (DR-RF)...\n")
tau.hat <- stage_2_rf(X, nuisances_rf$po_matrix, fold_indices, fold_list)

# step 2: calibrate sample.fraction
cat("Calibrating sample.fraction...\n")
cal <- find_optimal_sf(
  X            = X,
  Y            = Y,
  W            = W,
  nuisances_rf = nuisances_rf,
  tau.hat      = tau.hat,
  fold_indices = fold_indices,
  sf_grid      = seq(0.05, 0.5, 0.05),
  n_sim        = 50,
  CI_boot      = 100,
  alpha        = alpha
)

# step 3: final CIs using the calibrated sf
cat("Running bootstrap CIs with optimal sf =", cal$optimal_sf, "...\n")
final_ci <- rf_half_boot(
  X            = X,
  Y            = Y,
  W            = W,
  po           = nuisances_rf$po_matrix,
  tau          = tau.hat,
  CI_boot      = CI_boot,
  CI_sf        = cal$optimal_sf,
  alpha        = alpha,
  fold_indices = fold_indices,
  fold_list    = fold_list
)

warnings()

# save results
results <- list(
  tau            = tau.hat,
  hb_lb          = final_ci$hb_lb,
  hb_ub          = final_ci$hb_ub,
  optimal_sf     = cal$optimal_sf,
  coverage_curve = cal$coverage_curve,
  truth          = gen$truth,
  data           = data
)

output_dir <- file.path(dirname(path), "results/confidence_intervals/continuous/sf_calibration",
                         paste0("scenario_", scenario), n)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(output_dir, paste0("res_sim_", run, ".RDS")))

print("Simulation completed!")
