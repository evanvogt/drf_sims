##########
# title: extra - adding in the HTE tests to the missing data scenarios
##########

# libraries


# paths
path <- here()

# libraries
library(here)
library(dplyr)

# paths
path <- here()

# functions
source(here("missing/continuous/cts_miss_dgms.R"))
source(here("missing/continuous/cts_miss_models.R"))
source(here("utils.R"))

# simulation parameters
i <- as.numeric(commandArgs(trailingOnly = T))

n_folds <- 10
workers <- 2

params <- expand.grid(
  scenario = c(1:5),
  n = c(500),
  type = c("prognostic", "predictive", "both"),
  prop = c(0.3),
  mechanism = c("MAR", "AUX", "AUX-Y"), # only for 1 scenario?
  method = c("complete_cases", "mean_imputation", "missforest", "regression", 
             "missing_indicator", "IPW", "multiple_imputation", "none"),
  run = c(1:100),
  stringsAsFactors = F)

# redundant scenarios
params <- params %>%
  filter(!(scenario == 1 & (type != "prognostic" | mechanism == "AUX-Y")))

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

# read in existing result file
res_path <- file.path(dirname(here()),"results", "missing", "continuous",
                      paste0("scenario_", scenario), n, type, prop, mechanism,
                      method, paste0("res_sim_", run, ".RDS"))

result <- readRDS(res_path)

models <- c("causal_forest", "dr_random_forest", "dr_oracle", "dr_semi_oracle", "dr_superlearner")

models_run <- intersect(models, names(result))

# parameters needed for the testing routines

data <- result$data
Y <- data$Y
W <- data$W
X <- data[,-c(1,2)]

W_05 <- rep(0.5, nrow(X))
nuisances_rf <- result$nuisances_rf
nuisances_sl <- result$nuisances_sl

# NOTE: nuisances_rf$W.hat / Y0.hat / po are whole-sample OOB vectors, not
# double-crossfit matrices, as of the "oob + oob(s)" strategy change to
# R/cate_models.R::nuisance_rf - this script only works against result files
# saved after that change. Older saved RDS files carried *_matrix fields and
# needed rowMeans(..., na.rm = T) here instead.
for (model in models_run) {
  if (model == "causal_forest") {
    result$causal_forest$BLP_whole <- run_blp_whole(Y, W, nuisances_rf$W.hat,
                                                    nuisances_rf$Y0.hat, result$causal_forest$tau)

    result$causal_forest$independence_cate <- run_independence_test_whole(X, result$causal_forest$tau)
    result$causal_forest$independence_po <- run_independence_test_whole(X, nuisances_rf$po)
  }
  if (model == "dr_random_forest") {
    result$dr_random_forest$BLP_whole <- run_blp_whole(Y, W, nuisances_rf$W.hat,
                                                    nuisances_rf$Y0.hat, result$dr_random_forest$tau)

    result$dr_random_forest$independence_cate <- run_independence_test_whole(X, result$dr_random_forest$tau)
    result$dr_random_forest$independence_po <- run_independence_test_whole(X, nuisances_rf$po)
  }
  if (model == "dr_oracle") {
    result$dr_oracle$BLP_whole <- run_blp_whole(Y, W, W05, result$dr_oracle$Y0.hat, result$dr_oracle$tau)

    result$dr_oracle$independence_cate <- run_independence_test_whole(X, result$dr_oracle$tau)
    result$dr_oracle$independence_po <- run_independence_test_whole(X, nuisances_rf$po)

  }
}