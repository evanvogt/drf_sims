##########
# title: half-sample bootstrap estimation - correlated binary outcome
##########
# The retired confidence_intervals/binary/bin_ci_analysis.R (tag
# independent-ss-final) on the correlated sets: the same estimators and
# bootstrap (bin_corr_ci_models.R, formerly
# confidence_intervals/binary/bin_ci_models.R), with rho read from the grid row.

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
source(here("sample_size", "correlated", "confidence_intervals", "binary",
            "bin_corr_ci_config.R"))

# simulation parameters
i <- as.numeric(commandArgs(trailingOnly = T))

CI_boot <- 200
alpha <- 0.05
n_folds <- 10
workers <- 2

# The parameter grid lives in the study config, so this script and the
# check/collect scripts cannot disagree about what index i means.
param <- study$grid[i, ]
print(param)

scenario <- param$scenario
n <- param$n
CI_sf <- param$CI_sf
run <- param$run
rho <- param$rho

# set up simulation seed - the run alone, so each rho (and CI_sf) sees the same draws
setup_rng_stream(run)

# data generation
gen <- generate_binary_scenario_data(scenario, n, rho)

data <- gen$dataset

fmla_info <- get_binary_oracle_info(scenario, gen$bW, rho)

# fixed covariate-grid query points, as the parent study - see
# R/dgm_scenarios.R::build_query_grid. The grid does not depend on rho (the
# truth function tau(x) is the parent's), but at rho = 0.5 its corners are
# low-density points: see ../README.md, "The query grid".
set <- corr_set("binary", rho)
Z_query <- build_query_grid(scenario, set = set,
                            covariate_names = names(data)[-c(1, 2)])
grid_truth <- build_query_grid_truth(scenario, set = set, gen$bW, Z_query)


# Set up parallelisation
metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)


# run the methods and get CIs
results <- run_all_cate_methods(
  data = data,
  n_folds = n_folds,
  fmla_info = fmla_info,
  CI_boot = CI_boot,
  CI_sf = CI_sf,
  alpha = alpha,
  Z_query = Z_query
)
warnings()

results$data <- data
results$truth <- gen$truth
results$Z_query <- Z_query
results$grid_truth <- grid_truth

# Save results, under the path check/collect expect (combo_dir())
output_dir <- combo_dir(study, param)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(output_dir, paste0("res_sim_", run, ".RDS")))

print("Simulation completed!")
