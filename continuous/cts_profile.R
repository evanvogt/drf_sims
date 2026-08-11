##########
# title: resource profiling for the continuous-outcome array job
##########
# Runs one complete replicate through the same orchestrator cts_analysis.R uses
# (run_all_cate_methods -> R/cate_models.R's cate_methods), under a sweep of
# (n, workers, grf_threads), and records timings, memory and CPU so the PBS
# directives on cts_1.sh can be measured rather than guessed. Parallel to
# crossfitting/cf_profile.R; see that file's header for the syrup/multisession
# rationale, not repeated here.
#
#   Rscript continuous/cts_profile.R 1
#
# The specific thing this settles: cts_1.sh has always requested
# ompthreads=2 alongside cts_analysis.R's hardcoded workers <- 2, but
# num.threads was never set anywhere in R/cate_models.R - so grf defaulted to
# all visible cores per multisession worker regardless of what was allocated.
# cate_methods() and cts_analysis.R now accept/forward num.threads; this sweep
# measures whether the 2x2 pairing was actually oversubscribed once that knob
# does something, the same question cf_profile.R answered for crossfitting.
#
# Unlike crossfitting (fixed n = 500), continuous's single array job spans
# n in {100, 250, 500, 1000} under one #PBS -l line, so n is swept here too -
# sizing must cover the n = 1000 cell, not just a middle value.

library(dplyr)
library(furrr)
library(grf)
library(GenericML)
library(SuperLearner)
library(here)

if (!requireNamespace("syrup", quietly = TRUE)) {
  stop("the syrup package is required for profiling - install it into sim-env ",
       "with install.packages('syrup')", call. = FALSE)
}

# path
path <- here()

# functions
source(here("R", "utils.R"))           # setup_rng_stream, timed
source(here("continuous", "cts_dgms.R"))
source(here("continuous", "cts_models.R"))  # run_all_cate_methods -> cate_methods

# profiling parameters
i <- as.numeric(commandArgs(trailingOnly = T))

sample_interval <- 1   # seconds between syrup snapshots

prof_params <- expand.grid(
  scenario = c(1, 7),      # fewest / most covariates (no X3-X5 vs X3+X4+X5)
  n = c(100, 1000),        # smallest / largest n in the study grid
  run = c(1, 2),           # replicate to replicate variability
  workers = c(1, 2, 4),
  grf_threads = c(1, 2),
  stringsAsFactors = F
)

param <- prof_params[i,]
print(param)

scenario <- param$scenario
n <- param$n
run <- param$run
workers <- param$workers
grf_threads <- param$grf_threads

# same n-dependent rules as cts_analysis.R, so profiled cost matches production
n_folds <- dplyr::case_when(n == 100 ~ 4L, n == 250 ~ 5L, TRUE ~ 10L)
sl_lib <- if (n <= 100) {
  c("SL.glm", "SL.glmnet", "SL.gam", "SL.mean")
} else {
  c("SL.glm", "SL.glmnet", "SL.earth", "SL.gam", "SL.mean", "SL.ranger")
}

# set up simulation seed
setup_rng_stream(run)

gen <- generate_continuous_scenario_data(scenario, n)
fmla_info <- get_continuous_oracle_info(scenario, gen$bW)

# multisession workers are new R processes and inherit this, so setting it here
# does control their OpenMP thread pools even though this process's own libraries
# have already initialised - matches cf_profile.R / cts_analysis.R
Sys.setenv(OMP_NUM_THREADS = grf_threads)

metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)

parent_pid <- Sys.getpid()

# syrup() returns the measurement tibble, not the value of the expression - but it
# evaluates the expression in this environment, so `total` and `res` persist.
# verbose_timing = TRUE asks cate_methods() for a per-model timing breakdown
# (results$timings), used below to build arm_times.
usage <- syrup::syrup({
  total <- system.time({
    res <- run_all_cate_methods(
      data = gen$dataset,
      n_folds = n_folds,
      sl_lib = sl_lib,
      fmla_info = fmla_info,
      num.threads = grf_threads,
      verbose_timing = TRUE
    )
  })
}, interval = sample_interval)

plan(sequential)

# syrup reports every R process `ps` can see, which on a shared node can include
# another of your own jobs - keep only this tree, as cf_profile.R does
usage <- usage %>%
  filter(pid == parent_pid | ppid == parent_pid) %>%
  mutate(role = if_else(pid == parent_pid, "parent", "worker"),
         rss = as.numeric(rss),
         vms = as.numeric(vms))

# future falls back to sequential when workers == 1
n_procs <- n_distinct(usage$pid)
expected_procs <- if (workers > 1) workers + 1L else 1L
if (n_procs != expected_procs) {
  warning(sprintf("expected %d processes in the tree but saw %d",
                  expected_procs, n_procs))
}

by_snapshot <- usage %>%
  group_by(id, time) %>%
  summarise(total_rss = sum(rss, na.rm = TRUE),
            total_pct_cpu = sum(pct_cpu, na.rm = TRUE),
            .groups = "drop")

peak_rss_gb <- max(by_snapshot$total_rss, na.rm = TRUE) / 1024^3

by_process <- usage %>%
  group_by(pid, role) %>%
  summarise(max_rss_gb = max(rss, na.rm = TRUE) / 1024^3,
            median_pct_cpu = median(pct_cpu, na.rm = TRUE),
            max_pct_cpu = max(pct_cpu, na.rm = TRUE),
            .groups = "drop")

allocated_pct_cpu <- workers * grf_threads * 100
busy <- by_snapshot$total_pct_cpu[by_snapshot$total_pct_cpu > 10]
median_pct_cpu <- if (length(busy) > 0) median(busy) else NA_real_

# per-model timings - one row per top-level block cate_methods() timed (shared
# nuisance_rf overhead plus each estimator). Unlike crossfitting's arms, these
# don't split into nuisance/stage2 sub-times, since cate_methods doesn't expose
# that split per model.
arm_times <- tibble::tibble(arm = names(res$timings),
                            time_total = unlist(res$timings, use.names = FALSE))

profile <- list(
  index = i,
  scenario = scenario,
  n = n,
  run = run,
  workers = workers,
  grf_threads = grf_threads,
  n_folds = n_folds,
  n_models = length(res$timings),
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
  expected_procs = expected_procs,
  peak_rss_gb = peak_rss_gb,
  median_pct_cpu = median_pct_cpu,
  allocated_pct_cpu = allocated_pct_cpu,
  sample_interval = sample_interval,
  warnings = names(warnings()),
  available_cores = as.integer(parallelly::availableCores()),
  omp_num_threads = Sys.getenv("OMP_NUM_THREADS"),
  r_version = R.version.string,
  timestamp = Sys.time()
)

cat("\n--- profile summary ---\n")
cat(sprintf("n=%d workers=%d grf_threads=%d  elapsed=%.1fs  cpu-seconds=%.1f\n",
            n, workers, grf_threads, profile$elapsed, profile$cpu_seconds))
cat(sprintf("processes tracked: %d (expected %d%s)\n", n_procs, expected_procs,
            if (workers == 1) ", future runs sequentially at workers=1" else ""))
cat(sprintf("peak memory: %.2f gb summed across the tree\n", peak_rss_gb))
cat(sprintf("cpu while busy: %.0f%% observed against %.0f%% allocated%s\n",
            median_pct_cpu, allocated_pct_cpu,
            if (!is.na(median_pct_cpu) && median_pct_cpu > allocated_pct_cpu * 1.1) {
              "  <-- OVERSUBSCRIBED"
            } else ""))
cat("\nper process:\n")
print(as.data.frame(by_process), digits = 3, row.names = FALSE)
cat("\ntime by model (nuisance_rf is shared overhead, not a model of its own):\n")
print(as.data.frame(arrange(arm_times, desc(time_total))), digits = 3, row.names = FALSE)

if (length(profile$warnings) > 0) {
  cat(sprintf("\n%d distinct warnings, most frequent:\n", length(profile$warnings)))
  print(head(sort(table(profile$warnings), decreasing = TRUE), 5))
}

output_dir <- file.path(dirname(path), "results", "continuous", "profiling")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(profile, file.path(output_dir, paste0("prof_", i, ".RDS")))

print("Profiling run completed!")
