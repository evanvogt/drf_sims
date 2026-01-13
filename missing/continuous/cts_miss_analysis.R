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
source(here("utils.R"))

# parameters for the simulation
i <- as.numeric(commandArgs(trailingOnly = T))

workers <- 1

params <- expand.grid(scenario = c(1:5),
                      n = c(500),
                      prop = c(0.3),
                      type = c("prognostic", "predictive", "both"),
                      mechanism = c("MAR", "AUX"),# , "AUX-Y"), # only for 1 scenario
                      method = c("complete_cases", "mean_imputation", "missforest", "regression", "missing_indicator", "IPW", "none"),
                      run = c(1:100))

# can't have predictive missingness in the scenario where there are no predictive variables
params <- params %>%
  filter(!(scenario == 1 & type != "prognostic"))

# select parameters for this run
param <- params[i,]
scenario <- as.numeric(param$scenario)
n <- as.numeric(param$n)
prop <- as.numeric(param$prop)
type <- as.character(param$type)
mechanism <- as.character(param$mechanism)
method <- as.character(param$method)
run <- as.numeric(param$run)

sl_lib <- c("SL.glm", "SL.glmnet", "SL.earth", "SL.gam", "SL.mean", "SL.randomForest")

# set up simulation seed
setup_rng_stream(run)

# data generation and missing data handling

gen <- generate_and_process_continuous_data(scenario, n, TRUE, type, prop, mechanism, method)

data <- gen$dataset

n_folds <- ifelse(nrow(data) < 201, 5, 10)

fmla_info <- get_continuous_oracle_info(scenario, gen$bW)


# Run all the CATE models

results <- run_all_cate_methods(data, n_folds, workers, sl_lib, fmla_info, ipw = if (method == "IPW") gen$ipw else NULL)
warnings()

results$data <- data
results$truth <- gen$truth

if (param$method == "IPW") {
  results$ipw <- gen$ipw
}

# Save results
output_dir <- file.path("/rds/general/project/nihr_drf_simulations", "live", "results", "missing", "continuous", paste0("scenario_", scenario), n, type, prop, mechanism, method)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, here(output_dir, paste0("res_sim_", run, ".RDS")))

print("Simulation completed!")

