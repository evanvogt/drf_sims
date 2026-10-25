##########
# Title: competing-risk SNR check - DR-learners vs pseudo-value CF
##########
# competing_risk/surv_analysis.R on the grid row's scenario, with the X3 effect
# scaled by hte_mult and only the arms in SNR_MODELS (snr_config.R). Same arguments:
#   Rscript snr_analysis.R <i> [workers] [grf_threads]

# Report failures on STDOUT, with the call stack - see surv_analysis.R for why.
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
library(here)

# R/cate_models.R before surv_models.R, as in surv_analysis.R, so the parent
# study's definitions win where both exist
source(here("R", "utils.R"))
source(here("R", "cate_models.R"))
source(here("competing_risk", "surv_dgm.R"))
source(here("competing_risk", "surv_models.R"))
source(here("competing_risk", "snr", "snr_config.R"))

# Simulation parameters
args <- as.numeric(commandArgs(trailingOnly = T))
i <- args[1]
# workers/grf_threads default to 2/1, as surv_analysis.R
workers <- if (length(args) >= 2 && !is.na(args[2])) args[2] else 2
grf_threads <- if (length(args) >= 3 && !is.na(args[3])) args[3] else 1

horizon <- survival_scenario_params$event_horizon[1]

param <- study$grid[i, ]
print(param)

scenario <- param$scenario
n <- param$n
hte_mult <- param$hte_mult
run <- param$run

# only the scf arms split, and none is run here; kept so the fold_indices saved
# with the results match surv_analysis.R's
n_folds <- ifelse(n < 300, 5, 10)
t0 <- Sys.time()
# the run alone, as the main study: the n = 500 cells of a run are paired
setup_rng_stream(run)

# keeps rfsrc's main-process OOB fits and the BLAS in step with grf's
# num.threads; set before plan() so the workers start with it
Sys.setenv(OMP_NUM_THREADS = grf_threads)

metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)

gen <- generate_surv_data(
  scenario = scenario,
  n = n,
  rho = SNR_RHO,
  censoring = SNR_CENSORING,
  hte_mult = hte_mult
)

data <- gen$dataset

results <- all_cate_surv_models(
  data = data,
  n_folds = n_folds,
  horizon = horizon,
  num.threads = grf_threads,
  models = SNR_MODELS,
  estimands = SNR_ESTIMANDS
)
t1 <- Sys.time()
results$data <- data
results$truth <- gen$truth
print(t1 - t0)

# save results, under the path check/collect expect (combo_dir())
output_dir <- combo_dir(study, param)
save_run(results, output_dir, run)

print("Simulation completed!")
