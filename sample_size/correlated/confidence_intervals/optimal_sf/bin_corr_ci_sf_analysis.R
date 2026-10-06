##########
# title: optimal sample.fraction calibration simulation - correlated binary outcome
##########
# The retired confidence_intervals/optimal_sf/bin_ci_sf_analysis.R (tag
# independent-ss-final) on the correlated sets, with the grid read from the
# study config and rho from the grid row. The calibration settings are the
# parent's. find_optimal_sf() is in R/bootstrap_ci.R.

# libraries
library(here)
library(dplyr)

# path
path <- here()

# functions
source(here("R", "utils.R"))
source(here("sample_size", "correlated", "binary", "bin_corr_dgms.R"))
source(here("sample_size", "correlated", "confidence_intervals", "binary",
            "bin_corr_ci_models.R"))
source(here("R", "bootstrap_ci.R"))
source(here("sample_size", "correlated", "confidence_intervals", "optimal_sf",
            "bin_corr_ci_sf_config.R"))

# simulation parameters
args <- as.numeric(commandArgs(trailingOnly = T))
i <- args[1]

CI_boot <- 200
alpha   <- 0.1
workers <- if (length(args) >= 2 && !is.na(args[2])) args[2] else 2

# no CI_sf axis - that is what we are finding
param    <- study$grid[i, ]
print(param)

scenario <- param$scenario
n        <- param$n
run      <- param$run
rho      <- param$rho

# set up simulation seed - the run alone, so each rho sees the same draws
setup_rng_stream(run)

# data generation
gen  <- generate_binary_scenario_data(scenario, n, rho)
data <- gen$dataset

X     <- as.matrix(data[, -c(1:2)])
Y     <- data$Y
W     <- data$W
n_obs <- nrow(X)

# fixed covariate-grid query points, for the grid-based (as opposed to
# per-unit) band on the final CI below - see R/dgm_scenarios.R::build_query_grid
# and ../README.md, "The query grid".
set         <- corr_set("binary", rho)
Z_query     <- build_query_grid(scenario, set = set,
                                covariate_names = names(data)[-c(1, 2)])
grid_truth  <- build_query_grid_truth(scenario, set = set, gen$bW, Z_query)
Z_query_mat <- as.matrix(Z_query)

# set up parallelisation
metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)

# step 1: lightweight fit - nuisances + tau.hat only, no bootstrap CIs
cat("Computing nuisance functions...\n")
nuisances_rf <- nuisance_rf(X, Y, W)

cat("Estimating tau (DR-RF)...\n")
s1           <- stage2_whole_rf(X, nuisances_rf$po, Z_query = Z_query_mat)
tau.hat      <- s1$tau
tau_grid.hat <- s1$tau_grid

# step 2: calibrate sample.fraction
cat("Calibrating sample.fraction...\n")
cal <- find_optimal_sf(
  X            = X,
  Y            = Y,
  W            = W,
  nuisances_rf = nuisances_rf,
  tau.hat      = tau.hat,
  sf_grid      = seq(0.05, 0.5, 0.05),
  n_sim        = 50,
  CI_boot      = 100,
  alpha        = alpha
)

# step 3: final CIs using the calibrated sf
cat("Running bootstrap CIs with optimal sf =", cal$optimal_sf, "...\n")
final_ci <- rf_oob_half_boot(
  X        = X,
  Y        = Y,
  W        = W,
  po       = nuisances_rf$po,
  tau      = tau.hat,
  CI_boot  = CI_boot,
  CI_sf    = cal$optimal_sf,
  alpha    = alpha,
  Z_query  = Z_query_mat,
  tau_grid = tau_grid.hat
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
  data           = data,
  tau_grid       = tau_grid.hat,
  grid_lb        = final_ci$grid_lb,
  grid_ub        = final_ci$grid_ub,
  Z_query        = Z_query,
  grid_truth     = grid_truth
)

# under the path check/collect expect (combo_dir())
output_dir <- combo_dir(study, param)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(output_dir, paste0("res_sim_", run, ".RDS")))

print("Simulation completed!")
