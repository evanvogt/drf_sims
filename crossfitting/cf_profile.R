##########
# title: resource profiling for the crossfitting array job
##########
# Runs one complete replicate through the same orchestrator cf_analysis.R uses,
# under a sweep of (workers, grf_threads), and records timings and memory so the
# PBS directives on cf_1.sh can be measured rather than guessed.
#
# The specific thing being settled: continuous/jobscripts/cts_1.sh asks for
# ompthreads=2 while cts_analysis.R sets workers <- 2. grf is OpenMP-threaded and
# defaults to num.threads = NULL (all cores), so each multisession worker can claim
# 2 threads against 2 allocated cores. num.threads is not set anywhere in the repo.

library(dplyr)
library(furrr)
library(grf)
library(SuperLearner)
library(here)

# path
path <- here()

# functions
source(here("crossfitting", "cf_models.R"))

# profiling parameters
i <- as.numeric(commandArgs(trailingOnly = T))

n <- 500
n_folds <- 10L
n_test <- 2000

prof_params <- expand.grid(
  scenario = c(1, 6),      # fewest / most covariates (7 vs 9 columns)
  run = c(1:3),            # replicate to replicate variability
  workers = c(1, 2, 4),
  grf_threads = c(1, 2),
  stringsAsFactors = F
)

param <- prof_params[i,]
print(param)

scenario <- param$scenario
run <- param$run
workers <- param$workers
grf_threads <- param$grf_threads

sl_lib <- c("SL.glm", "SL.glmnet", "SL.earth", "SL.gam", "SL.mean", "SL.ranger")

# set up simulation seed
setup_rng_stream(run)

gen <- generate_cf_replicate(scenario, n, n_test)

# R level peak memory. this is a lower bound: multisession workers are separate
# processes and are invisible to this process's gc(). the true figure comes from
# the PBS accounting sampled by cf_profile.sh.
gc(reset = TRUE, full = TRUE)

# multisession workers are new R processes and inherit this, so setting it here
# does control their OpenMP thread pools even though this process's own libraries
# have already initialised. keeps SL.ranger and the BLAS in step with grf's
# num.threads instead of each claiming every core.
Sys.setenv(OMP_NUM_THREADS = grf_threads)

metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)

total <- system.time({
  res <- run_all_crossfit_variants(
    data = gen$data,
    X_test = gen$X_test,
    n_folds = n_folds,
    sl_lib = sl_lib,
    num.threads = grf_threads,
    truth_test = gen$truth_test_tau
  )
})

gc_after <- gc(full = TRUE)
plan(sequential)

# per arm timings
arm_times <- bind_rows(lapply(names(res$arms), function(nm) {
  a <- res$arms[[nm]]
  tibble::tibble(arm = nm, family = a$family, variant = a$variant,
                 time_nuisance = a$time_nuisance, time_stage2 = a$time_stage2,
                 time_total = a$time_nuisance + a$time_stage2)
}))

profile <- list(
  index = i,
  scenario = scenario,
  run = run,
  workers = workers,
  grf_threads = grf_threads,
  n = n,
  n_folds = n_folds,
  n_arms = length(res$arms),
  elapsed = unname(total["elapsed"]),
  user = unname(total["user.self"]),
  system = unname(total["sys.self"]),
  # cpu-seconds charged: what PBS bills, and what the throttled array is limited by
  cpu_seconds = unname(total["elapsed"]) * workers * grf_threads,
  arm_times = arm_times,
  # gc() columns are used/(Mb)/gc trigger/(Mb)/max used/(Mb) - the "(Mb)" names
  # are duplicated, so column 6 (peak, in Mb) is taken by position
  r_peak_mb = sum(gc_after[, 6]),
  gc_table = gc_after,
  available_cores = parallelly::availableCores(),
  omp_num_threads = Sys.getenv("OMP_NUM_THREADS"),
  r_version = R.version.string,
  timestamp = Sys.time()
)

cat("\n--- profile summary ---\n")
cat(sprintf("workers=%d grf_threads=%d  elapsed=%.1fs  cpu-seconds=%.1f\n",
            workers, grf_threads, profile$elapsed, profile$cpu_seconds))
print(arrange(arm_times, desc(time_total)))

output_dir <- file.path(dirname(path), "results", "crossfitting", "profiling")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(profile, file.path(output_dir, paste0("prof_", i, ".RDS")))

print("Profiling run completed!")
