##########
# title: resource profiling for the crossfitting array job
##########
# Runs one complete replicate through the same orchestrator cf_analysis.R uses,
# under a sweep of (workers, grf_threads), and records timings, memory and CPU so
# the PBS directives on cf_1.sh can be measured rather than guessed.
#
# Instrumented with syrup, which spins up a separate R session that samples `ps`
# at a fixed interval and reports every R process - so it sees the multisession
# workers, which the parent's own gc() cannot. Nothing here depends on the
# scheduler, so this script also runs on a laptop:
#
#   Rscript crossfitting/cf_profile.R 1
#
# The specific thing being settled: continuous/jobscripts/cts_1.sh asks for
# ompthreads=2 while cts_analysis.R sets workers <- 2. grf is OpenMP-threaded and
# defaults to num.threads = NULL (all cores), so each multisession worker can claim
# 2 threads against 2 allocated cores. num.threads is not set anywhere in the repo.
# syrup's pct_cpu measures that directly rather than inferring it from wall time.

library(dplyr)
library(furrr)
library(grf)
library(SuperLearner)
library(here)

if (!requireNamespace("syrup", quietly = TRUE)) {
  stop("the syrup package is required for profiling - install it into sim-env ",
       "with install.packages('syrup')", call. = FALSE)
}

# path
path <- here()

# functions
source(here("crossfitting", "cf_models.R"))

# profiling parameters
i <- as.numeric(commandArgs(trailingOnly = T))

n <- 500
n_folds <- 10L
n_test <- 2000
sample_interval <- 1   # seconds between syrup snapshots

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

# multisession workers are new R processes and inherit this, so setting it here
# does control their OpenMP thread pools even though this process's own libraries
# have already initialised. keeps SL.ranger and the BLAS in step with grf's
# num.threads instead of each claiming every core.
Sys.setenv(OMP_NUM_THREADS = grf_threads)

metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)

parent_pid <- Sys.getpid()

# syrup() returns the measurement tibble, not the value of the expression - but it
# evaluates the expression in this environment, so `total` and `res` persist
usage <- syrup::syrup({
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
}, interval = sample_interval)

plan(sequential)

# syrup reports every R process `ps` can see, which on a shared node can include
# another of your own jobs. multisession workers are PSOCK children of this
# process, so keep only this tree.
usage <- usage %>%
  filter(pid == parent_pid | ppid == parent_pid) %>%
  mutate(role = if_else(pid == parent_pid, "parent", "worker"))

n_procs <- n_distinct(usage$pid)
if (n_procs != workers + 1) {
  warning(sprintf("expected %d processes in the tree (1 parent + %d workers) but saw %d",
                  workers + 1, workers, n_procs))
}

# peak memory: the single worst snapshot, summed across the tree. this overcounts
# shared library pages, which is the safe direction for a mem= request.
by_snapshot <- usage %>%
  group_by(id, time) %>%
  summarise(total_rss = sum(rss, na.rm = TRUE),
            total_pct_cpu = sum(pct_cpu, na.rm = TRUE),
            .groups = "drop")

peak_rss_gb <- max(by_snapshot$total_rss, na.rm = TRUE) / 1024^3

# per process summary - shows whether the workers are balanced and the parent idle
by_process <- usage %>%
  group_by(pid, role) %>%
  summarise(max_rss_gb = max(rss, na.rm = TRUE) / 1024^3,
            median_pct_cpu = median(pct_cpu, na.rm = TRUE),
            max_pct_cpu = max(pct_cpu, na.rm = TRUE),
            .groups = "drop")

# observed against allocated: this is the oversubscription measurement
allocated_pct_cpu <- workers * grf_threads * 100
# only count snapshots where work was actually happening
busy <- by_snapshot$total_pct_cpu[by_snapshot$total_pct_cpu > 10]
median_pct_cpu <- if (length(busy) > 0) median(busy) else NA_real_

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
  usage = usage,                 # full time series, for the timeline plot
  by_snapshot = by_snapshot,
  by_process = by_process,
  n_procs = n_procs,
  peak_rss_gb = peak_rss_gb,
  median_pct_cpu = median_pct_cpu,
  allocated_pct_cpu = allocated_pct_cpu,
  sample_interval = sample_interval,
  available_cores = as.integer(parallelly::availableCores()),
  omp_num_threads = Sys.getenv("OMP_NUM_THREADS"),
  r_version = R.version.string,
  timestamp = Sys.time()
)

cat("\n--- profile summary ---\n")
cat(sprintf("workers=%d grf_threads=%d  elapsed=%.1fs  cpu-seconds=%.1f\n",
            workers, grf_threads, profile$elapsed, profile$cpu_seconds))
cat(sprintf("processes tracked: %d (expected %d)\n", n_procs, workers + 1))
cat(sprintf("peak memory: %.2f gb summed across the tree\n", peak_rss_gb))
cat(sprintf("cpu while busy: %.0f%% observed against %.0f%% allocated%s\n",
            median_pct_cpu, allocated_pct_cpu,
            if (!is.na(median_pct_cpu) && median_pct_cpu > allocated_pct_cpu * 1.1) {
              "  <-- OVERSUBSCRIBED"
            } else ""))
cat("\nper process:\n")
print(as.data.frame(by_process), digits = 3, row.names = FALSE)
cat("\nslowest arms:\n")
print(as.data.frame(head(arrange(arm_times, desc(time_total)), 8)), digits = 3, row.names = FALSE)

output_dir <- file.path(dirname(path), "results", "crossfitting", "profiling")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(profile, file.path(output_dir, paste0("prof_", i, ".RDS")))

print("Profiling run completed!")
