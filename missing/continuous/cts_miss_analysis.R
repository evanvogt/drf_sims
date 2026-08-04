##########
# title: fit all CATE models under various missingness settings and handlings
##########

# libraries
library(here)
library(dplyr)

# paths
path <- here()

# functions
source(here("missing/continuous/cts_miss_dgms.R"))
source(here("missing/continuous/cts_miss_models.R"))
source(here("missing/continuous/cts_miss_config.R"))
source(here("utils.R"))

# simulation parameters
i <- as.numeric(commandArgs(trailingOnly = T))

n_folds <- 10
workers <- 2

# The grid lives in cts_miss_config.R and is NEVER filtered here. It used to be
# (`filter(method == "complete_data")`), which renumbered every row, so index i
# meant one thing here and another in cts_miss_check.R - a failed_ids.txt would
# have resubmitted the wrong parameters. To run one arm, select indices instead:
#   idx <- grid_indices(study, method = "complete_data")
# and submit those, leaving the numbering intact.
param <- study$grid[i, ]
print(param)

scenario <- param$scenario
n <- param$n
prop <- param$prop
type <- param$type
mechanism <- param$mechanism
method <- param$method
run <- param$run

sl_lib <- c("SL.glm", "SL.glmnet", "SL.earth", "SL.gam", "SL.mean", "SL.ranger")

# set up simulation seed
setup_rng_stream(run)

if (method == "complete_data") {
  # Reference run: complete data without any missingness, for relative efficiency benchmarking
  gen <- generate_continuous_scenario_data(scenario, n, mech = mechanism, return_truth = TRUE)
  data <- gen$dataset
  fmla_info <- get_continuous_oracle_info(scenario, gen$bW)

  metaplan <- plan(multisession, workers = workers)
  on.exit(plan(metaplan), add = TRUE)

  results <- run_all_cate_methods(
    data = data,
    n_folds = n_folds,
    sl_lib = sl_lib,
    fmla_info = fmla_info
  )
  warnings()
} else {
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

  # set up parallelisation
  metaplan <- plan(multisession, workers = workers)
  on.exit(plan(metaplan), add = TRUE)

  # Run all the CATE models
  if (method == "multiple_imputation") {
    result_list <- future_map(data, function(data) {
      run_all_cate_methods(
        data = data,
        n_folds = n_folds,
        sl_lib = NULL,
        fmla_info = NULL
      )
    }, .options = furrr_options(seed = TRUE))

    warnings()
    results <- list()
    results$causal_forest <- combine_mi(result_list, "causal_forest")
    results$dr_random_forest <- combine_mi(result_list, "dr_random_forest")
    results$dr_semi_oracle <- combine_mi(result_list, "dr_semi_oracle")
  } else {
    results <- run_all_cate_methods(
      data = data,
      n_folds = n_folds,
      sl_lib = sl_lib,
      fmla_info = fmla_info,
      ipw = if (method == "IPW") gen$ipw else NULL
    )
    warnings()
  }
  plan(metaplan)
}

results$data <- data
results$truth <- gen$truth

if (param$method == "IPW") {
  results$ipw <- gen$ipw
}

# Save results
out_path <- file.path(dirname(here()),"results", "missing", "continuous", paste0("scenario_", scenario), n, type, prop, mechanism, method)
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(out_path, paste0("res_sim_", run, ".RDS")))

print("Simulation completed!")
