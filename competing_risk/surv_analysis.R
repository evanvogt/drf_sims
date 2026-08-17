##########
# Title: Competing risks CATEs
##########

# Report failures on STDOUT, not stderr.
#
# The cluster is currently keeping only the PBS `.o` files, and R sends errors -
# and every message() in all_cate_surv_models() - to stderr, so when 225 array
# indices died there was nothing in the logs to say why and it took a local
# reproduction to find out (see surv_failed_diagnose.R). This handler puts the
# condition message and the call stack where PBS will keep them. It still exits
# non-zero, so a failed job is still a failed job as far as the scheduler and
# check_failed() are concerned.
#
# It is armed before the library()/source() calls, and uses only base functions,
# so a missing package or an unparseable source file is reported the same way.
options(error = function() {
  calls <- sys.calls()
  calls <- calls[-length(calls)] # drop this handler's own frame
  cat("\n=== SIMULATION FAILED ===\n")
  cat("args:", paste(commandArgs(trailingOnly = TRUE), collapse = " "), "\n")
  cat("error:", trimws(geterrmessage()), "\n")
  cat("call stack (outermost first):\n")
  for (k in seq_along(calls)) {
    cat(sprintf("  %2d: %s\n", k, deparse(calls[[k]])[1]))
  }
  cat("=========================\n")
  flush(stdout())
  quit(status = 1, save = "no")
})

library(dplyr)
library(furrr)
library(grf)
library(GenericML)
library(SuperLearner)
library(here)

# Path
path <- here()

# Functions
# R/cate_models.R is sourced BEFORE surv_models.R so that this study's own
# definitions win where they still exist. It supplies the shared crossfitting
# machinery this study now shares with the rest of the repo:
# oob_predict_counterfactual, stage2_whole_rf, stage_2_sl, pretest_superlearner
# and dr_pseudo. trim_ps arrives via R/utils.R, which cate_models.R sources.
source(here("R", "utils.R"))
source(here("R", "cate_models.R"))
source(here("competing_risk", "surv_dgm.R"))
source(here("competing_risk", "surv_models.R"))
source(here("competing_risk/surv_config.R"))

# Simulation parameters
i <- as.numeric(commandArgs(trailingOnly = T))

workers <- 2

horizon <- 28

# The parameter grid lives in the study config, so this script and the
# check/collect scripts cannot disagree about what index i means.
param <- study$grid[i, ]
print(param)

scenario <- param$scenario
n <- param$n
censoring <- param$censoring
run <- param$run

n_folds <- ifelse(n < 300, 5, 10)
t0 <- Sys.time()
# Set up simulation seed
setup_rng_stream(run)

# Dataset Generation
metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)

gen <- generate_surv_data(
  scenario = scenario,
  n = n,
  censoring = censoring
)

data <- gen$dataset

# Main analysis
results <- all_cate_surv_models(
  data = data,
  n_folds = n_folds,
  horizon = horizon
)
t1 <- Sys.time()
results$data <- data
results$truth <- gen$truth
print(t1-t0)
# save results
output_dir <- file.path(dirname(path), "results", "competing_risk", paste0("scenario_", scenario), n, paste0("censor_", censoring))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(output_dir, paste0("res_sim_", run, ".RDS")))

print("Simulation completed!")