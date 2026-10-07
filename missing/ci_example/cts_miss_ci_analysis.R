##########
# title: confidence intervals example for multiple imputation + bootstrap
##########

# libraries
library(here)

# path
path <- here()

# functions
source(here("missing", "ci_example", "cts_miss_ci_dgms.R"))
source(here("missing", "ci_example", "cts_miss_ci_models.R"))
source(here("R", "utils.R"))
source(here("missing/ci_example/cts_miss_ci_config.R"))

# simulation parameters. Trailing args after the array index are
# workers/grf_threads, as the jobscripts pass them. workers is the number of
# imputations fitted in parallel; everything inside an imputation (crossfit
# folds, bootstrap draws) runs sequentially. Without them the defaults are
# 3 workers and grf's own thread default.
args <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(args[1])

n_folds <- 10
workers <- if (length(args) >= 2) as.numeric(args[2]) else 3
grf_threads <- if (length(args) >= 3) as.numeric(args[3]) else NULL
CI_boot <- 200
CI_sf <- 0.5
alpha <- 0.05

# The parameter grid lives in the study config, so this script and the
# check/collect scripts cannot disagree about what index i means.
param <- study$grid[i, ]
print(param)

scenario <- param$scenario
n <- param$n
prop <- param$prop
type <- param$type
mechanism <- param$mechanism
method <- param$method
run <- param$run

# set up simulation seed
setup_rng_stream(run)

# data generation and missing data handling
gen <- generate_and_process_continuous_data(
  scenario = scenario,
  n = n, 
  return_truth = TRUE, 
  type = type,
  prop = prop, 
  mech = mechanism,
  method = method)

data <- gen$dataset

fmla_info <- get_continuous_oracle_info(scenario, gen$bW)

# set up parallelisation. multisession workers are new R processes and
# inherit this, so setting it here controls their OpenMP thread pools even
# though this process's own libraries have already initialised.
# Only the imputation level is parallel: a nested multisession level made each
# worker open its own sockets, which failed on the cluster when array jobs
# sharing a node picked the same port. With one worker nothing opens a socket.
if (!is.null(grf_threads)) Sys.setenv(OMP_NUM_THREADS = grf_threads)
outer_plan <- if (workers > 1) tweak(multisession, workers = workers) else sequential
metaplan <- plan(list(outer_plan, sequential))
on.exit(plan(metaplan), add = T)

# Run the models
results <- mi_boot(
  datalist = data,
  n_folds = n_folds,
  fmla_info = fmla_info,
  CI_boot = CI_boot,
  CI_sf = CI_sf,
  alpha = alpha,
  num.threads = grf_threads
)
plan(metaplan)
results$data <- data
results$truth <- gen$truth
# which covariates the amputation left missing (n x 5), for coverage on the
# complete and incomplete units separately (cts_miss_ci_metrics.R)
results$miss_mask <- gen$miss_mask
results$retained_indices <- gen$retained_indices

# Save results
out_path <- file.path(dirname(here()), "results", "missing", "ci_example", paste0("scenario_", scenario), n, type, prop, mechanism, method)
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(out_path, paste0("res_sim_", run, ".RDS")))

print("Simulation completed!")