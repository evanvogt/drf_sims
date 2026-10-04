##########
# Title: single-event survival CATEs
##########

# Report failures on STDOUT, with the call stack - see competing_risk/
# surv_analysis.R, which this copies. It exits non-zero, so a failed job is
# still a failed job as far as the scheduler and check_failed() are concerned.
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
library(SuperLearner)
library(here)

# Functions
# R/cate_models.R is sourced BEFORE competing_risk/surv_models.R, as in
# surv_analysis.R, so the parent study's definitions win where both exist.
# se_models.R reuses the parent's arms and adds the single-event pieces.
source(here("R", "utils.R"))
source(here("R", "cate_models.R"))
source(here("competing_risk", "surv_models.R"))
source(here("competing_risk", "single_event", "se_dgm.R"))
source(here("competing_risk", "single_event", "se_models.R"))
source(here("competing_risk", "single_event", "se_config.R"))

# Simulation parameters
args <- as.numeric(commandArgs(trailingOnly = T))
i <- args[1]
# workers/grf_threads default to 2/1, as surv_analysis.R
workers <- if (length(args) >= 2 && !is.na(args[2])) args[2] else 2
grf_threads <- if (length(args) >= 3 && !is.na(args[3])) args[3] else 1

horizon <- se_scenario_params$event_horizon[1]

param <- study$grid[i, ]
print(param)

scenario <- param$scenario
n <- param$n
censoring <- param$censoring
run <- param$run

n_folds <- ifelse(n < 300, 5, 10)
t0 <- Sys.time()
# Set up simulation seed - the run alone, so the scenarios and censoring
# settings of a run share their draws (and W, X and U match competing_risk's
# rho = 0 run of the same index - see se_dgm.R)
setup_rng_stream(run)

# keeps rfsrc's main-process OOB fits and the BLAS in step with grf's
# num.threads; set before plan() so the workers start with it
Sys.setenv(OMP_NUM_THREADS = grf_threads)

metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)

gen <- generate_se_data(
  scenario = scenario,
  n = n,
  censoring = censoring
)

data <- gen$dataset

# Main analysis
results <- all_cate_se_models(
  data = data,
  n_folds = n_folds,
  horizon = horizon,
  sl_library = sl_libraries(n),
  num.threads = grf_threads
)
t1 <- Sys.time()
results$data <- data
results$truth <- gen$truth
print(t1 - t0)

# save results, under the path check/collect expect (combo_dir())
output_dir <- combo_dir(study, param)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(output_dir, paste0("res_sim_", run, ".RDS")))

print("Simulation completed!")
