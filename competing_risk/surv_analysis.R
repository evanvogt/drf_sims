##########
# Title: Competing risks CATEs
##########

# Report failures on STDOUT, with the call stack.
#
# Added when the cluster was keeping only the PBS `.o` files: R sends errors -
# and every message() in all_cate_surv_models() - to stderr, so when 225 array
# indices died there was nothing in the logs to say why and it took a local
# reproduction to find out (see surv_failed_diagnose.R). The jobscripts now use
# `#PBS -j oe`, so stderr is merged into the `.o` file and R's own "Error in"
# line lands there too (the message appears twice). The handler is kept for
# what R does not print by default: the array args and the call stack. It still
# exits non-zero, so a failed job is still a failed job as far as the scheduler
# and check_failed() are concerned.
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
# t_learner_rf, stage2_whole_rf, stage_2_sl, pretest_superlearner
# and dr_pseudo. trim_ps arrives via R/utils.R, which cate_models.R sources.
source(here("R", "utils.R"))
source(here("R", "cate_models.R"))
source(here("competing_risk", "surv_dgm.R"))
source(here("competing_risk", "surv_models.R"))
source(here("competing_risk/surv_config.R"))

# Simulation parameters
args <- as.numeric(commandArgs(trailingOnly = T))
i <- args[1]
# workers/grf_threads default to 2/1, so `Rscript surv_analysis.R <i>` - what
# surv_1.sh, surv_2.sh and surv_run.R run - keeps 2 workers and now pins grf to
# 1 thread, the same arg order and defaults as sample_size/correlated/continuous/cts_corr_analysis.R.
workers <- if (length(args) >= 2 && !is.na(args[2])) args[2] else 2
grf_threads <- if (length(args) >= 3 && !is.na(args[3])) args[3] else 1

horizon <- 28

# The parameter grid lives in the study config, so this script and the
# check/collect scripts cannot disagree about what index i means.
param <- study$grid[i, ]
print(param)

scenario <- param$scenario
n <- param$n
censoring <- param$censoring
run <- param$run
rho <- param$rho

n_folds <- ifelse(n < 300, 5, 10)
t0 <- Sys.time()
# Set up simulation seed - the run alone, so each rho sees the same draws
setup_rng_stream(run)

# multisession workers are new R processes and inherit this, so it keeps rfsrc's
# main-process OOB fits and the BLAS in step with grf's num.threads rather than
# each process claiming every core - matches sample_size/correlated/continuous/cts_corr_analysis.R.
# Set before plan() so the workers start with it.
Sys.setenv(OMP_NUM_THREADS = grf_threads)

# Dataset Generation
metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)

gen <- generate_surv_data(
  scenario = scenario,
  n = n,
  rho = rho,
  censoring = censoring
)

data <- gen$dataset

# Main analysis
results <- all_cate_surv_models(
  data = data,
  n_folds = n_folds,
  horizon = horizon,
  sl_library = sl_libraries(n),
  num.threads = grf_threads
)
t1 <- Sys.time()
results$data <- data
results$truth <- gen$truth
print(t1-t0)
# save results, under the path check/collect expect (combo_dir())
output_dir <- combo_dir(study, param)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(output_dir, paste0("res_sim_", run, ".RDS")))

print("Simulation completed!")