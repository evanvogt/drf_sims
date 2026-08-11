##########
# title: resource profiling for the missing-binary array job
##########
# Runs one complete replicate through the same orchestrator bin_miss_analysis.R
# uses (generate_and_process_binary_data / generate_binary_scenario_data ->
# run_all_cate_methods -> R/cate_models.R's cate_methods), under a sweep of
# (method, scenario, mechanism, run, workers, grf_threads), and records
# timings, memory and CPU so the PBS directives on bin_miss_1.sh can be
# measured rather than guessed. Parallel to binary/bin_profile.R and
# continuous/cts_profile.R; see those files' headers for the syrup/multisession
# rationale, not repeated here.
#
#   Rscript missing/binary/bin_miss_profile.R 1
#
# The specific thing this settles: bin_miss_1.sh has always requested
# ompthreads=2 alongside bin_miss_analysis.R's hardcoded workers <- 2, but
# num.threads was never set anywhere in bin_miss_models.R - so grf defaulted
# to all visible cores per multisession worker regardless of what was
# allocated. cate_methods(), bin_miss_models.R and bin_miss_analysis.R now
# accept/forward num.threads; this sweep measures whether that pairing was
# actually oversubscribed once the knob does something, the same question
# bin_profile.R answered for binary/.
#
# TWO differences from bin_profile.R/cts_profile.R, both driven by this
# study's method axis (bin_miss_config.R's MISS_METHODS_STUDY, 9 values):
#
# 1. Data generation is NOT free here. generate_and_process_binary_data()
#    calls mice()/missForest()/etc *before* plan(multisession) is even set up
#    - it is serial, single-threaded, and independent of (workers,
#    grf_threads). Re-running it once per (workers, grf_threads) cell (as a
#    flat expand.grid sweep would) would burn most of the profiling budget
#    re-measuring the same number six times over. Instead this script uses a
#    two-level grid: an OUTER grid of (method, scenario, mechanism, run) - one
#    PBS array index each - generates/processes data ONCE per outer row
#    (timed with plain system.time(), no syrup needed for a single serial
#    process), then an INNER grid of (workers, grf_threads) loops over that
#    same processed data, syrup-wrapping only the parallel model-fitting
#    phase per inner cell, exactly as bin_profile.R does for its one cell.
#    Peak memory is not measured separately for the data-gen phase: the R
#    process persists across both phases within one task, so peak_rss_gb
#    measured during the (first) inner cell already reflects whatever
#    mice/missForest output is still retained in scope - no double counting.
#
# 2. method == "multiple_imputation" doesn't call run_all_cate_methods() once
#    - bin_miss_analysis.R's own multiple_imputation branch calls future_map()
#    over the 50 datasets mice() produced, then combine_mi()s three of the
#    five arms. That whole block (future_map + combine_mi) is what gets
#    syrup-wrapped for that method's inner cells, not a single
#    run_all_cate_methods() call. combine_mi() does not carry $timings
#    through, so this script sums each arm's per-imputation elapsed time
#    across the 50 runs into one arm_times table, comparable in shape (if not
#    in what it represents) to the other methods' single-call arm_times.

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
source(here("R", "utils.R"))                       # setup_rng_stream, timed
source(here("missing", "binary", "bin_miss_dgms.R"))
source(here("missing", "binary", "bin_miss_models.R"))  # run_all_cate_methods -> cate_methods, combine_mi
source(here("missing", "binary", "bin_miss_config.R"))  # MISS_METHODS_STUDY

# profiling parameters
i <- as.numeric(commandArgs(trailingOnly = T))

sample_interval <- 1   # seconds between syrup snapshots

# outer grid: one PBS array index each - the expensive, workers/grf_threads-
# independent data generation + missingness handling happens once per row.
# scenario 1/4 are the fewest/most-covariates extremes within {1,2,4,5} (the
# study's actual scenario set), per R/dgm_scenarios.R's binary_missing table
# (scenario 1: no X3/X4/X5; scenario 4: X3+X4+X5). mechanism drops MNAR-Y -
# mechanistically closest to MNAR for a resource sweep, and excluded for
# scenario 1 in bin_miss_config.R's own grid anyway.
outer_params <- expand.grid(
  method = MISS_METHODS_STUDY,
  scenario = c(1, 4),
  mechanism = c("MAR", "MNAR"),
  run = c(1, 2),
  stringsAsFactors = FALSE
)

# inner grid: looped inside this one task, reusing the outer row's data
inner_params <- expand.grid(
  workers = c(1, 2, 4),
  grf_threads = c(1, 2),
  stringsAsFactors = FALSE
)

outer <- outer_params[i, ]
print(outer)

scenario <- outer$scenario
mechanism <- outer$mechanism
method <- outer$method
run <- outer$run

# fixed by bin_miss_config.R's grid - this study does not sweep n/type/prop
n <- 500
type <- "both"
prop <- 0.3
n_folds <- 10
sl_lib <- c("SL.glm", "SL.glmnet", "SL.earth", "SL.gam", "SL.mean", "SL.ranger")

# set up simulation seed - same stream a production run at this (scenario,
# mechanism, method, run) combination would use
setup_rng_stream(run)

# ---- data generation + missingness handling: serial, one-off ---------------
data_gen_time <- system.time({
  if (method == "complete_data") {
    gen <- generate_binary_scenario_data(scenario, n, mech = mechanism, return_truth = TRUE)
  } else {
    gen <- generate_and_process_binary_data(
      scenario = scenario,
      n = n,
      return_truth = TRUE,
      type = type,
      prop = prop,
      mech = mechanism,
      method = method)
  }
})

fmla_info <- get_binary_oracle_info(scenario, gen$bW)

cat(sprintf("\ndata generation + missingness handling (method=%s): %.1fs\n",
            method, data_gen_time["elapsed"]))

parent_pid <- Sys.getpid()

# ---- inner sweep: parallel model-fitting phase, reusing the same data ------
inner_results <- vector("list", nrow(inner_params))

for (j in seq_len(nrow(inner_params))) {

  workers <- inner_params$workers[j]
  grf_threads <- inner_params$grf_threads[j]
  cat(sprintf("\n--- inner cell %d/%d: workers=%d grf_threads=%d ---\n",
              j, nrow(inner_params), workers, grf_threads))

  # multisession workers are new R processes and inherit this, so setting it
  # here does control their OpenMP thread pools even though this process's
  # own libraries have already initialised - matches bin_profile.R
  Sys.setenv(OMP_NUM_THREADS = grf_threads)

  metaplan <- plan(multisession, workers = workers)

  # syrup::syrup() is called up to 6 times per array task here, unlike
  # bin_profile.R/cts_profile.R which call it exactly once - and on this local
  # Windows machine (see memory: syrup profiling is unreliable here) that
  # repeated use has been observed to fail deterministically on a later call
  # in the loop ("Unable to retrieve resource usage results from the temporary
  # session", a syrup-internal error) after several earlier calls in the same
  # session succeeded. Wrapping each cell in tryCatch means one bad
  # measurement costs that cell, not the other 5 - the array task still
  # completes and saves whatever it collected, rather than losing an entire
  # (method, scenario, mechanism, run) row's worth of profiling for one
  # syrup hiccup. plan(sequential) is reset in `finally` so a mid-cell failure
  # can't leave a stale multisession plan for the next iteration.
  cell <- tryCatch({

    # syrup() returns the measurement tibble, not the value of the expression -
    # but it evaluates the expression in this environment, so `total` and `res`
    # persist. Mirrors bin_profile.R's single-cell wrap; the branching inside
    # matches bin_miss_analysis.R's three method branches exactly.
    usage <- syrup::syrup({
      total <- system.time({
        if (method == "complete_data") {
          res <- run_all_cate_methods(
            data = gen$dataset,
            n_folds = n_folds,
            sl_lib = sl_lib,
            fmla_info = fmla_info,
            num.threads = grf_threads,
            verbose_timing = TRUE
          )
        } else if (method == "multiple_imputation") {
          result_list <- future_map(gen$dataset, function(d) {
            run_all_cate_methods(
              data = d,
              n_folds = n_folds,
              sl_lib = NULL,
              fmla_info = NULL,
              num.threads = grf_threads,
              verbose_timing = TRUE
            )
          }, .options = furrr_options(seed = TRUE))

          res <- list()
          res$causal_forest <- combine_mi(result_list, "causal_forest")
          res$dr_random_forest <- combine_mi(result_list, "dr_random_forest")
          res$dr_semi_oracle <- combine_mi(result_list, "dr_semi_oracle")

          # combine_mi() drops $timings - sum each arm's elapsed time across the
          # 50 imputations so this method's arm_times is comparable in shape to
          # the other methods', even though it represents a 50x sum, not a
          # single call
          timing_names <- unique(unlist(lapply(result_list, function(x) names(x$timings))))
          res$timings <- setNames(
            lapply(timing_names, function(nm) {
              sum(vapply(result_list, function(x) {
                v <- x$timings[[nm]]
                if (is.null(v)) 0 else v
              }, numeric(1)))
            }),
            timing_names
          )
        } else {
          res <- run_all_cate_methods(
            data = gen$dataset,
            n_folds = n_folds,
            sl_lib = sl_lib,
            fmla_info = fmla_info,
            ipw = if (method == "IPW") gen$ipw else NULL,
            num.threads = grf_threads,
            verbose_timing = TRUE
          )
        }
      })
    }, interval = sample_interval)

    # syrup reports every R process `ps` can see, which on a shared node can
    # include another of your own jobs - keep only this tree, as bin_profile.R does
    usage <- usage %>%
      filter(pid == parent_pid | ppid == parent_pid) %>%
      mutate(role = if_else(pid == parent_pid, "parent", "worker"),
             rss = as.numeric(rss),
             vms = as.numeric(vms))

    # future falls back to sequential when workers == 1
    n_procs <- n_distinct(usage$pid)
    expected_procs <- if (workers > 1) workers + 1L else 1L
    if (n_procs != expected_procs) {
      warning(sprintf("cell %d: expected %d processes in the tree but saw %d",
                      j, expected_procs, n_procs))
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

    arm_times <- tibble::tibble(arm = names(res$timings),
                                time_total = unlist(res$timings, use.names = FALSE))

    list(
      ok = TRUE,
      workers = workers,
      grf_threads = grf_threads,
      n_models = length(res$timings),
      elapsed = unname(total["elapsed"]),
      user = unname(total["user.self"]),
      system = unname(total["sys.self"]),
      # cpu-seconds charged for the model-fitting phase only - data_gen_elapsed
      # (outer, above) is serial and single-threaded regardless of this cell
      cpu_seconds = unname(total["elapsed"]) * workers * grf_threads,
      arm_times = arm_times,
      usage = usage,
      by_snapshot = by_snapshot,
      by_process = by_process,
      n_procs = n_procs,
      expected_procs = expected_procs,
      peak_rss_gb = peak_rss_gb,
      median_pct_cpu = median_pct_cpu,
      allocated_pct_cpu = allocated_pct_cpu
    )
  }, error = function(e) {
    warning(sprintf("cell %d (workers=%d, grf_threads=%d) failed: %s",
                    j, workers, grf_threads, conditionMessage(e)))
    list(ok = FALSE, workers = workers, grf_threads = grf_threads,
         error = conditionMessage(e))
  }, finally = {
    plan(sequential)
  })

  inner_results[[j]] <- cell
}

profile <- list(
  index = i,
  method = method,
  scenario = scenario,
  mechanism = mechanism,
  run = run,
  n = n,
  n_folds = n_folds,
  data_gen_elapsed = unname(data_gen_time["elapsed"]),
  data_gen_user = unname(data_gen_time["user.self"]),
  inner = inner_results,
  sample_interval = sample_interval,
  warnings = names(warnings()),
  available_cores = as.integer(parallelly::availableCores()),
  r_version = R.version.string,
  timestamp = Sys.time()
)

cat("\n--- profile summary ---\n")
cat(sprintf("method=%s scenario=%d mechanism=%s run=%d\n", method, scenario, mechanism, run))
cat(sprintf("data generation + handling: %.1fs\n", profile$data_gen_elapsed))
cat("\ninner cells (parallel model-fitting phase only; NA rows are cells where syrup itself failed):\n")
inner_df <- purrr::map_dfr(inner_results, function(x) {
  data.frame(workers = x$workers, grf_threads = x$grf_threads, ok = x$ok,
             elapsed = if (isTRUE(x$ok)) x$elapsed else NA_real_,
             cpu_seconds = if (isTRUE(x$ok)) x$cpu_seconds else NA_real_,
             peak_rss_gb = if (isTRUE(x$ok)) x$peak_rss_gb else NA_real_,
             n_procs = if (isTRUE(x$ok)) x$n_procs else NA_integer_,
             expected_procs = if (isTRUE(x$ok)) x$expected_procs else NA_integer_,
             median_pct_cpu = if (isTRUE(x$ok)) x$median_pct_cpu else NA_real_,
             allocated_pct_cpu = if (isTRUE(x$ok)) x$allocated_pct_cpu else NA_real_)
})
print(inner_df, digits = 3, row.names = FALSE)
if (all(!inner_df$ok)) {
  warning("every inner cell failed for this outer row - see the per-cell warnings above")
} else {
  cat(sprintf("\ntotal replicate walltime at each config = data_gen_elapsed + elapsed, e.g. best case %.1fs, worst case %.1fs\n",
              profile$data_gen_elapsed + min(inner_df$elapsed, na.rm = TRUE),
              profile$data_gen_elapsed + max(inner_df$elapsed, na.rm = TRUE)))
}

if (length(profile$warnings) > 0) {
  cat(sprintf("\n%d distinct warnings, most frequent:\n", length(profile$warnings)))
  print(head(sort(table(profile$warnings), decreasing = TRUE), 5))
}

output_dir <- file.path(dirname(path), "results", "missing", "binary", "profiling")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(profile, file.path(output_dir, paste0("prof_", i, ".RDS")))

print("Profiling run completed!")
