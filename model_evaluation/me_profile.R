##########
# title: resource profiling for the model evaluation array job
##########
# Runs one complete replicate through the same two-phase orchestration
# me_analysis.R uses (candidate models, then nuisance-evaluation pipelines),
# under a sweep of (workers, n_cores), and records timings, memory and CPU so
# the PBS directives on me_1.sh can be measured rather than guessed - the
# same mechanism crossfitting/cf_profile.R uses for cf_1.sh.
#
# Instrumented with syrup, which spins up a separate R session that samples
# `ps` at a fixed interval and reports every R process - so it sees the
# multisession workers, which the parent's own gc() cannot. Nothing here
# depends on the scheduler, so this also runs on a laptop as a smoke test:
#
#   Rscript model_evaluation/me_profile.R 1
#
# Deliberately smaller sweep than cf_profile.R's 36 cells (2 scenarios x 3
# runs x 3 workers x 2 thread-values): each cell here is a real H2O AutoML +
# XGBoost-CV run, not a handful of RF fits, so this sweeps 2 scenarios x 2
# runs x 2 worker-values x 2 n_cores-values = 16 cells, at the smallest grid
# n (250) to keep each cell as fast as possible. Neither knob includes a
# value of 1 in this sweep (unlike cf_profile.R's grf_threads/workers), so
# there is no workers==1/n_cores==1 baseline to compute a speedup-efficiency
# table against - me_profile_summary.R reports absolute timings instead.
#
# Candidate-model fitting (future-parallel, controlled by `workers`) and
# nuisance evaluation (XGBoost/H2O-parallel, controlled by `n_cores`) run
# SEQUENTIALLY in me_analysis.R, governed by different knobs - so they are
# timed separately here (unlike cf_profile.R, which times one combined
# call). A combined number would hide which phase dominates, which is
# exactly the information the resource decision needs.
#
# CAVEAT SPECIFIC TO THIS STUDY: syrup tracks R processes via ps (pid ==
# parent_pid | ppid == parent_pid). H2O's own JVM is a separate Java
# process, and whether it appears as a tracked child of this filter or is
# invisible to it is not yet confirmed. This script prints the full,
# unfiltered process list syrup observed (not just an expected-count
# assertion, the way cf_profile.R does) so that's visible in the first
# profiling run rather than silently under-counting H2O's memory footprint.

library(dplyr)
library(future)
library(future.apply)
library(ranger)
library(glmnet)
library(SuperLearner)
library(xgboost)
library(h2o)
library(caret)
library(here)

if (!requireNamespace("syrup", quietly = TRUE)) {
  stop("the syrup package is required for profiling - install it into sim-env ",
       "with install.packages('syrup')", call. = FALSE)
}

# path
path <- here()

# functions
source(here("R", "utils.R"))
source(here("model_evaluation", "me_dgms.R"))
source(here("model_evaluation", "me_utils.R"))
source(here("model_evaluation", "me_models.R"))
source(here("model_evaluation", "me_nuisance.R"))
source(here("model_evaluation", "me_config.R"))

# profiling parameters
i <- as.numeric(commandArgs(trailingOnly = T))

n <- 250 # smallest grid n
n_folds <- 10L # matches me_analysis.R's rule (single crossfitting, all n)
h2o_mem <- "10G"
sample_interval <- 1 # seconds between syrup snapshots

# SuperLearner is noisy - keep them all rather than losing everything past
# the 50th, and stash them in the profile object so they can be inspected
options(nwarnings = 500)

prof_params <- expand.grid(
  scenario = c(1, 6), # cheapest / priciest of this study's scenario set
  run = c(1, 2), # replicate-to-replicate variability
  workers = c(2, 4),
  n_cores = c(2, 5),
  stringsAsFactors = FALSE
)

param <- prof_params[i, ]
print(param)

scenario <- param$scenario
run <- param$run
workers <- param$workers
n_cores <- param$n_cores

# set up simulation seed
setup_rng_stream(run)

gen <- generate_me_scenario_data(scenario, n)
design <- prepare_design_matrix(gen$dataset)
Y <- design$Y
W <- design$W
X <- design$X
kfolds <- split_folds(Y, k = n_folds)
nuis_folds <- split_folds(Y, k = n_folds)

metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)

parent_pid <- Sys.getpid()

# syrup() returns the measurement tibble, not the value of the expression -
# but it evaluates the expression in this environment, so the system.time()
# results and the fitted objects persist afterward
usage <- syrup::syrup({
  total_models <- system.time({
    model_list <- run_all_candidate_models(
      Y, W, X, kfolds$fold_indices, kfolds$fold_list
    )
  })

  model_seed <- sample.int(2^31 - 1, 1)
  total_nuisances <- system.time({
    nuisances <- run_all_nuisance_pipelines(
      X, Y, W, nuis_folds$fold_indices, nuis_folds$fold_list,
      n_cores = n_cores, mem = h2o_mem, model_seed = model_seed
    )
  })
}, interval = sample_interval)

plan(sequential)

# the caveat this file's header describes: print every process syrup saw,
# before filtering, so an H2O JVM that isn't a tracked child of this R
# process is visible here rather than silently dropped by the filter below
cat("\n--- all processes syrup observed (before filtering to this job's R tree) ---\n")
print(as.data.frame(usage %>% distinct(pid, ppid)), row.names = FALSE)

# syrup reports every R process `ps` can see, which on a shared node can
# include another of your own jobs. multisession workers are PSOCK children
# of this process, so keep only this tree.
usage <- usage %>%
  filter(pid == parent_pid | ppid == parent_pid) %>%
  mutate(role = if_else(pid == parent_pid, "parent", "worker"),
         # syrup returns rss/vms as bench::bench_bytes so they print as
         # "1.2GB". that class survives arithmetic, so coerce here, before
         # anything divides.
         rss = as.numeric(rss),
         vms = as.numeric(vms))

# future falls back to sequential when workers == 1; not reachable in this
# sweep (workers is always >= 2 here), kept for consistency with cf_profile.R
n_procs <- n_distinct(usage$pid)
expected_procs <- if (workers > 1) workers + 1L else 1L
if (n_procs != expected_procs) {
  warning(sprintf(
    paste0("expected %d processes in the tree but saw %d after filtering to ",
           "pid/ppid == %d. If H2O's JVM is not a tracked child of this R ",
           "process, its memory is NOT included in peak_rss_gb below - see ",
           "the unfiltered process list printed above, and cross-check ",
           "against `qstat -fx <jobid> | grep resources_used` once this ",
           "runs for real."),
    expected_procs, n_procs, parent_pid
  ))
}

# peak memory: the single worst snapshot, summed across the tree. this
# overcounts shared library pages, which is the safe direction for a mem=
# request - modulo the H2O-visibility caveat above.
by_snapshot <- usage %>%
  group_by(id, time) %>%
  summarise(total_rss = sum(rss, na.rm = TRUE),
            total_pct_cpu = sum(pct_cpu, na.rm = TRUE),
            .groups = "drop")

peak_rss_gb <- max(by_snapshot$total_rss, na.rm = TRUE) / 1024^3

# per process summary - shows whether the workers are balanced and the
# parent idle
by_process <- usage %>%
  group_by(pid, role) %>%
  summarise(max_rss_gb = max(rss, na.rm = TRUE) / 1024^3,
            median_pct_cpu = median(pct_cpu, na.rm = TRUE),
            max_pct_cpu = max(pct_cpu, na.rm = TRUE),
            .groups = "drop")

# workers and n_cores parallelise two DIFFERENT, SEQUENTIAL phases (not
# concurrent), so there's no single "allocated" figure that's exactly right
# across the whole replicate - this is an approximation using whichever
# phase asks for more.
allocated_pct_cpu <- max(workers, n_cores) * 100
busy <- by_snapshot$total_pct_cpu[by_snapshot$total_pct_cpu > 10]
median_pct_cpu <- if (length(busy) > 0) median(busy) else NA_real_

profile <- list(
  index = i,
  scenario = scenario,
  run = run,
  workers = workers,
  n_cores = n_cores,
  n = n,
  n_folds = n_folds,
  elapsed_models = unname(total_models["elapsed"]),
  elapsed_nuisances = unname(total_nuisances["elapsed"]),
  elapsed = unname(total_models["elapsed"]) + unname(total_nuisances["elapsed"]),
  # cpu-seconds charged: what PBS bills, and what a throttled array is
  # limited by - each phase weighted by its OWN parallelism knob
  cpu_seconds = unname(total_models["elapsed"]) * workers +
    unname(total_nuisances["elapsed"]) * n_cores,
  usage = usage, # full time series, for the timeline plot
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
  r_version = R.version.string,
  timestamp = Sys.time()
)

cat("\n--- profile summary ---\n")
cat(sprintf(
  "workers=%d n_cores=%d  candidate-models=%.1fs  nuisances=%.1fs  total=%.1fs  cpu-seconds=%.1f\n",
  workers, n_cores, profile$elapsed_models, profile$elapsed_nuisances,
  profile$elapsed, profile$cpu_seconds
))
cat(sprintf("processes tracked: %d (expected %d)\n", n_procs, expected_procs))
cat(sprintf("peak memory: %.2f gb summed across the tracked R process tree\n", peak_rss_gb))
cat(sprintf("cpu while busy: %.0f%% observed against %.0f%% allocated (approx, see note above)%s\n",
            median_pct_cpu, allocated_pct_cpu,
            if (!is.na(median_pct_cpu) && median_pct_cpu > allocated_pct_cpu * 1.1) {
              "  <-- OVERSUBSCRIBED"
            } else ""))
cat("\nper process:\n")
print(as.data.frame(by_process), digits = 3, row.names = FALSE)

if (length(profile$warnings) > 0) {
  cat(sprintf("\n%d distinct warnings, most frequent:\n", length(profile$warnings)))
  print(head(sort(table(profile$warnings), decreasing = TRUE), 5))
}

output_dir <- file.path(dirname(path), "results", "model_evaluation", "profiling")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(profile, file.path(output_dir, paste0("prof_", i, ".RDS")))

print("Profiling run completed!")
