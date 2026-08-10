##########
# title: resource profiling for the CI pilot's half-sample bootstrap
##########
# Runs one complete cf_ci_analysis.R replicate - all 12 RF/CF arms plus their
# half-sample bootstrap CIs (R/bootstrap_ci.R's {rf,cf}_half_boot for the 5
# crossfit-structured arms, {rf,cf}_oob_half_boot for the 7 whole-sample/OOB
# ones) - under a sweep of (workers, grf_threads, CI_boot), and records
# timings, memory and CPU so cf_ci_1.sh's PBS directives can be measured rather
# than left at their placeholder values. Parallel to crossfitting/cf_profile.R,
# which does the same for the point-estimate study's cf_1.sh; see that file's
# header for the syrup/multisession rationale, not repeated here.
#
#   Rscript crossfitting/confidence_intervals/cf_ci_profile.R 1
#
# CI_boot is swept at two small values (20, 60) rather than at the pilot's real
# CI_boot = 200: each bootstrap draw is an independent future_map() iteration
# that refits V forests on a half-sample, and the n x CI_boot draws matrix is
# only assembled after future_map returns - so elapsed time should scale
# ~linearly in CI_boot while peak memory (governed by how many draws are being
# fit concurrently across `workers`, not by CI_boot) should be roughly flat.
# Profiling directly at CI_boot = 200 across this grid would cost as much as the
# production run it exists to size; two points let cf_ci_profile_summary.R fit a
# line, extrapolate to 200, and flag it if that linearity assumption doesn't
# actually hold.
#
# KNOWN WRONG, not fixed here: the "peak memory is roughly flat in CI_boot"
# claim above. future_map accumulates all CI_boot results before the draws matrix
# is assembled, so memory does grow with CI_boot - and cf_ci_profile_summary.R
# extrapolates only elapsed time, applying a flat mem_factor to peak RSS observed
# at CI_boot = 20/60. The mem= it writes into cf_ci_1.sh is therefore an
# underestimate of the pilot's real CI_boot = 200. See crossfitting/README.md.

library(dplyr)
library(furrr)
library(grf)
library(here)

if (!requireNamespace("syrup", quietly = TRUE)) {
  stop("the syrup package is required for profiling - install it into sim-env ",
       "with install.packages('syrup')", call. = FALSE)
}

# path
path <- here()

# functions
source(here("crossfitting", "cf_models.R"))   # generate_cf_replicate, run_all_crossfit_variants, timed
source(here("R", "bootstrap_ci.R"))           # {rf,cf}_half_boot, {rf,cf}_oob_half_boot

# profiling parameters
i <- as.numeric(commandArgs(trailingOnly = T))

n <- 500
n_folds <- 10L
n_test <- 2000
alpha <- 0.05
CI_sf <- 0.5             # cf_ci_analysis.R's fixed sample.fraction - no sweep here either
sample_interval <- 1     # seconds between syrup snapshots

prof_params <- expand.grid(
  scenario = c(1, 6),      # fewest / most covariates (7 vs 9 columns), same as cf_profile.R
  run = c(1, 2),           # replicate to replicate variability
  workers = c(1, 2, 4),
  grf_threads = c(1, 2),
  CI_boot = c(20, 60),     # two points to fit and sanity-check linear scaling in CI_boot
  stringsAsFactors = F
)

param <- prof_params[i,]
print(param)

scenario <- param$scenario
run <- param$run
workers <- param$workers
grf_threads <- param$grf_threads
CI_boot <- param$CI_boot

# set up simulation seed
setup_rng_stream(run)

gen <- generate_cf_replicate(scenario, n, n_test)

# multisession workers are new R processes and inherit this, so setting it here
# does control their OpenMP thread pools even though this process's own libraries
# have already initialised - matches cf_profile.R
Sys.setenv(OMP_NUM_THREADS = grf_threads)

metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)

parent_pid <- Sys.getpid()

# syrup() returns the measurement tibble, not the value of the expression - but it
# evaluates the expression in this environment, so `total`, `structured` and
# `boot_times` persist. The block covers the whole replicate cf_ci_analysis.R
# runs: nuisances + stage 2 for the 12 arms, then their bootstrap CIs in turn.
usage <- syrup::syrup({
  total <- system.time({
    structured <- run_all_crossfit_variants(
      data = gen$data,
      X_test = gen$X_test,
      n_folds = n_folds,
      sl_lib = NULL,
      num.threads = grf_threads,
      truth_test = gen$truth_test_tau
    )

    X <- as.matrix(gen$data[, -c(1:2)])
    Y <- gen$data$Y
    W <- gen$data$W
    nz <- structured$nuisances
    fold_indices <- structured$fold_indices
    fold_list <- unique(fold_indices)
    fold_indices_b <- structured$fold_indices_b
    fold_list_b <- unique(fold_indices_b)

    # table-driven bootstrap wiring, identical to cf_ci_analysis.R
    boot_spec <- list(
      # per-fold crossfit stage 2
      dcf            = list(fn = rf_half_boot,     arg = nz$nz_double$po,     fi = fold_indices,   fl = fold_list),
      scf_scf        = list(fn = rf_half_boot,     arg = nz$nz_single$po,     fi = fold_indices,   fl = fold_list),
      scf_scf_new    = list(fn = rf_half_boot,     arg = nz$nz_single$po,     fi = fold_indices_b, fl = fold_list_b),
      cf_dcf         = list(fn = cf_half_boot,     arg = nz$nz_double,        fi = fold_indices,   fl = fold_list),
      cf_scf         = list(fn = cf_half_boot,     arg = nz$nz_single,        fi = fold_indices,   fl = fold_list),
      # whole-sample stage 2, OOB predictions
      scf_oob        = list(fn = rf_oob_half_boot, arg = nz$nz_single$po,     fi = NULL,           fl = NULL),
      scf_oob_t      = list(fn = rf_oob_half_boot, arg = nz$nz_single_t$po,   fi = NULL,           fl = NULL),
      oob_oob        = list(fn = rf_oob_half_boot, arg = nz$nz_oob$po,        fi = NULL,           fl = NULL),
      oob_oob_s      = list(fn = rf_oob_half_boot, arg = nz$nz_oob_s$po,      fi = NULL,           fl = NULL),
      oob_oob_manual = list(fn = rf_oob_half_boot, arg = nz$nz_oob_manual$po, fi = NULL,           fl = NULL),
      cf_full_oob    = list(fn = cf_oob_half_boot, arg = nz$nz_single,        fi = NULL,           fl = NULL),
      cf_default     = list(fn = cf_oob_half_boot, arg = nz$nz_cf_default,    fi = NULL,           fl = NULL)
    )

    boot_times <- list()
    for (nm in names(boot_spec)) {
      spec <- boot_spec[[nm]]
      b <- timed(spec$fn(X, Y, W, spec$arg, structured$arms[[nm]]$tau,
                         CI_boot, CI_sf, alpha, spec$fi, spec$fl))
      boot_times[[nm]] <- b$time
      # draws intentionally dropped, as in cf_ci_analysis.R - only the timing matters here
    }
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

# nuisance/stage2 timings for all 12 arms - CI_boot-independent
arm_times <- bind_rows(lapply(names(structured$arms), function(nm) {
  a <- structured$arms[[nm]]
  tibble::tibble(arm = nm, family = a$family, variant = a$variant,
                 time_nuisance = a$time_nuisance, time_stage2 = a$time_stage2,
                 time_total = a$time_nuisance + a$time_stage2)
}))

# bootstrap timings for the same 12 arms - this is what scales with CI_boot
boot_times_tbl <- tibble::tibble(arm = names(boot_times),
                                 time_boot = unlist(boot_times, use.names = FALSE))

profile <- list(
  index = i,
  scenario = scenario,
  run = run,
  workers = workers,
  grf_threads = grf_threads,
  CI_boot = CI_boot,
  CI_sf = CI_sf,
  n = n,
  n_folds = n_folds,
  n_arms = length(structured$arms),
  elapsed = unname(total["elapsed"]),
  user = unname(total["user.self"]),
  system = unname(total["sys.self"]),
  # cpu-seconds charged: what PBS bills, and what the throttled array is limited by
  cpu_seconds = unname(total["elapsed"]) * workers * grf_threads,
  arm_times = arm_times,
  boot_times = boot_times_tbl,
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
cat(sprintf("workers=%d grf_threads=%d CI_boot=%d  elapsed=%.1fs  cpu-seconds=%.1f\n",
            workers, grf_threads, CI_boot, profile$elapsed, profile$cpu_seconds))
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
cat("\nslowest arms (nuisance + stage 2, CI_boot-independent):\n")
print(as.data.frame(head(arrange(arm_times, desc(time_total)), 8)), digits = 3, row.names = FALSE)
cat(sprintf("\nslowest bootstrap arms (at this cell's CI_boot=%d):\n", CI_boot))
print(as.data.frame(arrange(boot_times_tbl, desc(time_boot))), digits = 3, row.names = FALSE)

if (length(profile$warnings) > 0) {
  cat(sprintf("\n%d distinct warnings, most frequent:\n", length(profile$warnings)))
  print(head(sort(table(profile$warnings), decreasing = TRUE), 5))
}

output_dir <- file.path(dirname(path), "results", "crossfitting_ci", "profiling")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(profile, file.path(output_dir, paste0("prof_", i, ".RDS")))

print("Profiling run completed!")
