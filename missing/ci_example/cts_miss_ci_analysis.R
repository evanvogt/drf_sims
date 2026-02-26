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
source(here("utils.R"))

# simulation parameters
i <- as.numeric(commandArgs(trailingOnly = T))

n_folds <- 10
workers <- c(3, 3)
CI_boot <- 200
CI_sf <- 0.5
alpha <- 0.05

params <- expand.grid(
  scenario = c(1:5),
  n = c(500),
  type = c("both"),
  prop = c(0.3),
  mechanism = c("MAR"),
  method = c("multiple_imputation"),
  run = c(1:100),
  stringsAsFactors = F
)

params <- params %>%
  arrange(scenario, run)

# select parameters for this run
param <- params[i,]
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

# set up parallelisation
metaplan <- plan(list(tweak(multisession, workers = workers[1]), tweak(multisession, workers = I(workers[2]))))
on.exit(plan(metaplan), add = T)

# Run the models
results <- mi_boot(
  datalist = data,
  n_folds = n_folds,
  fmla_info = fmla_info,
  CI_boot = CI_boot,
  CI_sf = CI_sf,
  alpha = alpha
)
plan(metaplan)
results$data <- data
results$truth <- gen$truth

# Save results
out_path <- file.path(dirname(here()), "results", "missing", "ci_example", paste0("scenario_", scenario), n, type, prop, mechanism, method)
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(out_path, paste0("res_sim_", run, ".RDS")))

print("Simulation completed!")