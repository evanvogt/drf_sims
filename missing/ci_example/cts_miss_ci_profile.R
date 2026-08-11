##########
# title: resource profiling for the ci_example multiple-imputation + bootstrap array job
##########
# Runs one complete replicate through the same building blocks
# cts_miss_ci_analysis.R uses (generate_and_process_continuous_data ->
# run_all_cate_methods/combine_mi_ci, i.e. mi_boot()'s own logic reimplemented
# inline - see below for why), under a sweep of (workers1, workers2,
# grf_threads, CI_boot), and records timings, memory and CPU so
# jobscripts/cts_miss_ci.sh's PBS directives can be measured rather than
# guessed. Parallel to missing/continuous/cts_miss_profile.R and
# crossfitting/confidence_intervals/cf_ci_profile.R; see those files' headers
# for the syrup/multisession rationale, not repeated here.
#
#   Rscript missing/ci_example/cts_miss_ci_profile.R 1
#
# This study is the only one in the repo with a NESTED two-level future plan
# (cts_miss_ci_analysis.R's plan(list(tweak(multisession, workers=workers[1]),
# tweak(multisession, workers=I(workers[2]))))). Tracing mi_boot()
# (cts_miss_ci_models.R): level 1 (workers[1]) parallelises future_map() over
# the 50 multiply-imputed datasets; level 2 (workers[2], wrapped I() so future
# won't auto-downscale it) parallelises the CI_boot=200 half-sample bootstrap
# draws *inside* each of 4 CATE arms' cf_oob_half_boot()/rf_oob_half_boot()
# calls (R/bootstrap_ci.R), run sequentially one arm at a time per imputation.
# So a full replicate is up to 50 x 4 x 200 bootstrap-forest refits, gated by
# two independent worker counts - this sweep exists to find out whether the
# current workers=c(3,3)/ncpus=10 guess was ever right.
#
# TWO differences from cts_miss_profile.R/cf_ci_profile.R, both driven by this
# study's shape:
#
# 1. Outer/inner split, like cts_miss_profile.R, because
#    generate_and_process_continuous_data() (mice(m=50)) is serial and
#    independent of (workers1, workers2, grf_threads, CI_boot) - but UNLIKE
#    cts_miss_profile.R, method is fixed to "multiple_imputation" throughout
#    (cts_miss_ci_config.R's grid never varies it), so every single outer row
#    pays that expensive path; there is no cheap complete_data branch to
#    special-case. The outer grid here is just (scenario, run) - the study's
#    only two free grid columns.
#
# 2. Reduced-CI_boot-and-extrapolate, like cf_ci_profile.R, because profiling
#    at the real CI_boot=200 costs as much as the production run it exists to
#    size. CI_boot is crossed with every (workers1, workers2, grf_threads)
#    inner cell (not swept separately) so cts_miss_ci_profile_summary.R can
#    fit elapsed ~ CI_boot per configuration and extrapolate to 200.
#
# mi_boot()/run_all_cate_methods() now forward num.threads (added alongside
# this script) so grf_threads is controllable, but mi_boot() still pools
# results via combine_mi_ci() before returning - which drops $timings, same
# as combine_mi() does for the non-CI missing-data studies. So this script
# doesn't call mi_boot() itself; it reimplements its ~15-line body directly
# (future_map + combine_mi_ci), calling run_all_cate_methods(...,
# verbose_timing = TRUE) so each imputation's $timings survives long enough
# to sum across the 50 imputations - identical in spirit to how
# cts_miss_profile.R's own multiple_imputation branch bypasses its study's
# combine-and-pool flow for the same reason.

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
source(here("R", "utils.R"))                                # setup_rng_stream
source(here("missing", "ci_example", "cts_miss_ci_dgms.R"))
source(here("missing", "ci_example", "cts_miss_ci_models.R")) # run_all_cate_methods; combine_mi_ci
                                                               # reaches this script transitively via
                                                               # cate_models.R's own source() of bootstrap_ci.R

# profiling parameters
i <- as.numeric(commandArgs(trailingOnly = T))

sample_interval <- 1   # seconds between syrup snapshots

# outer grid: one PBS array index each - the expensive, (workers1, workers2,
# grf_threads, CI_boot)-independent data generation + mice(m=50) imputation
# happens once per row. All 5 scenarios (this study's whole range, already
# the "reduced set" per its own README) since data-gen is paid once per row,
# not multiplied by the inner grid - cheap to cover all of them and get a
# per-scenario cost breakdown out of it.
outer_params <- expand.grid(
  scenario = 1:5,
  run = c(1, 2),
  stringsAsFactors = FALSE
)

# inner grid: looped inside this one task, reusing the outer row's data.
# workers1/workers2 bracket the production default of c(3,3) with the
# sequential case; grf_threads and CI_boot are crossed with every
# (workers1, workers2) pair because the summary script fits elapsed ~
# CI_boot PER (workers1, workers2, grf_threads) configuration.
inner_params <- expand.grid(
  workers1 = c(1, 3),
  workers2 = c(1, 3),
  grf_threads = c(1, 2),
  CI_boot = c(20, 60),
  stringsAsFactors = FALSE
)

outer <- outer_params[i, ]
print(outer)

scenario <- outer$scenario
run <- outer$run

# fixed by cts_miss_ci_config.R's grid - this study does not sweep n/type/prop/mechanism/method
n <- 500
type <- "both"
prop <- 0.3
mechanism <- "MAR"
n_folds <- 10
alpha <- 0.05
CI_sf <- 0.5

# set up simulation seed - same stream a production run at this (scenario, run)
# combination would use
setup_rng_stream(run)

# ---- data generation + missingness handling: serial, one-off ---------------
data_gen_time <- system.time({
  gen <- generate_and_process_continuous_data(
    scenario = scenario,
    n = n,
    return_truth = TRUE,
    type = type,
    prop = prop,
    mech = mechanism,
    method = "multiple_imputation")
})

fmla_info <- get_continuous_oracle_info(scenario, gen$bW)

cat(sprintf("\ndata generation + mice(m=50) imputation (scenario=%d): %.1fs\n",
            scenario, data_gen_time["elapsed"]))

parent_pid <- Sys.getpid()

# ---- inner sweep: parallel model-fitting + bootstrap phase, reusing the same data ------
inner_results <- vector("list", nrow(inner_params))

for (j in seq_len(nrow(inner_params))) {

  workers1 <- inner_params$workers1[j]
  workers2 <- inner_params$workers2[j]
  grf_threads <- inner_params$grf_threads[j]
  CI_boot <- inner_params$CI_boot[j]
  cat(sprintf("\n--- inner cell %d/%d: workers1=%d workers2=%d grf_threads=%d CI_boot=%d ---\n",
              j, nrow(inner_params), workers1, workers2, grf_threads, CI_boot))

  # multisession workers are new R processes and inherit this, so setting it
  # here does control their OpenMP thread pools even though this process's
  # own libraries have already initialised - matches cts_miss_profile.R
  Sys.setenv(OMP_NUM_THREADS = grf_threads)

  plan(list(tweak(multisession, workers = workers1),
            tweak(multisession, workers = I(workers2))))

  # syrup::syrup() is called up to 16 times per array task here, unlike
  # bin_profile.R/cts_profile.R which call it exactly once - and on this local
  # Windows machine (see memory: syrup profiling is unreliable here) that
  # repeated use has been observed to fail deterministically on a later call
  # in the loop after several earlier calls in the same session succeeded.
  # Wrapping each cell in tryCatch means one bad measurement costs that cell,
  # not the other 15 - the array task still completes and saves whatever it
  # collected. plan(sequential) is reset in `finally` so a mid-cell failure
  # can't leave a stale nested plan for the next iteration.
  cell <- tryCatch({

    # syrup() returns the measurement tibble, not the value of the expression -
    # but it evaluates the expression in this environment, so `total` and
    # `res` persist.
    usage <- syrup::syrup({
      total <- system.time({
        # mi_boot()'s own body, reimplemented so per-imputation $timings
        # survive long enough to sum (mi_boot()/combine_mi_ci() drop them
        # during pooling - see this file's header)
        results_list <- future_map(gen$dataset, function(d) {
          run_all_cate_methods(
            data = d,
            n_folds = n_folds,
            fmla_info = fmla_info,
            CI_boot = CI_boot,
            CI_sf = CI_sf,
            alpha = alpha,
            num.threads = grf_threads,
            verbose_timing = TRUE
          )
        }, .options = furrr_options(seed = TRUE))

        res <- list()
        res$causal_forest <- combine_mi_ci(results_list, "causal_forest", alpha)
        res$dr_random_forest <- combine_mi_ci(results_list, "dr_random_forest", alpha)
        res$dr_oracle <- combine_mi_ci(results_list, "dr_oracle", alpha)
        res$dr_semi_oracle <- combine_mi_ci(results_list, "dr_semi_oracle", alpha)

        # combine_mi_ci() drops $timings - sum each arm's elapsed time across
        # the 50 imputations, identical pattern to cts_miss_profile.R's
        # multiple_imputation branch (there for combine_mi(), here for
        # combine_mi_ci())
        timing_names <- unique(unlist(lapply(results_list, function(x) names(x$timings))))
        res$timings <- setNames(
          lapply(timing_names, function(nm) {
            sum(vapply(results_list, function(x) {
              v <- x$timings[[nm]]
              if (is.null(v)) 0 else v
            }, numeric(1)))
          }),
          timing_names
        )
      })
    }, interval = sample_interval)

    # syrup reports every R process `ps` can see, which on a shared node can
    # include another of your own jobs - keep only this tree, as cts_miss_profile.R does
    usage <- usage %>%
      filter(pid == parent_pid | ppid == parent_pid) %>%
      mutate(role = if_else(pid == parent_pid, "parent", "worker"),
             rss = as.numeric(rss),
             vms = as.numeric(vms))

    # Hypothesis, not a certainty: plan(list(tweak(w1), tweak(I(w2)))) should
    # mean each level-1 worker independently realises its own level-2 pool
    # the first time it nests a future (cf_oob_half_boot()/rf_oob_half_boot()'s
    # own future_map()), so at peak: controller + workers1 level-1 workers +
    # workers1*workers2 level-2 workers. Flagged as a warning (not asserted)
    # precisely so this sweep's own n_procs can confirm or refute it, the
    # same spirit as the existing expected_procs checks in
    # cts_miss_profile.R/cf_ci_profile.R. future falls back to sequential
    # when workers==1 at either level, hence the three degenerate branches.
    n_procs <- n_distinct(usage$pid)
    expected_procs <- if (workers1 > 1 && workers2 > 1) {
      1L + workers1 + workers1 * workers2
    } else if (workers1 > 1) {
      1L + workers1
    } else if (workers2 > 1) {
      1L + workers2
    } else {
      1L
    }
    if (n_procs != expected_procs) {
      warning(sprintf(
        "cell %d: expected %d processes (hypothesis: 1 + workers1 + workers1*workers2) but saw %d",
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

    allocated_pct_cpu <- workers1 * workers2 * grf_threads * 100
    busy <- by_snapshot$total_pct_cpu[by_snapshot$total_pct_cpu > 10]
    median_pct_cpu <- if (length(busy) > 0) median(busy) else NA_real_

    arm_times <- tibble::tibble(arm = names(res$timings),
                                time_total = unlist(res$timings, use.names = FALSE))

    list(
      ok = TRUE,
      workers1 = workers1,
      workers2 = workers2,
      grf_threads = grf_threads,
      CI_boot = CI_boot,
      n_models = length(res$timings),
      elapsed = unname(total["elapsed"]),
      user = unname(total["user.self"]),
      system = unname(total["sys.self"]),
      # cpu-seconds charged for the model-fitting + bootstrap phase only -
      # data_gen_elapsed (outer, above) is serial and single-threaded
      # regardless of this cell
      cpu_seconds = unname(total["elapsed"]) * workers1 * workers2 * grf_threads,
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
    warning(sprintf("cell %d (workers1=%d, workers2=%d, grf_threads=%d, CI_boot=%d) failed: %s",
                    j, workers1, workers2, grf_threads, CI_boot, conditionMessage(e)))
    list(ok = FALSE, workers1 = workers1, workers2 = workers2, grf_threads = grf_threads,
         CI_boot = CI_boot, error = conditionMessage(e))
  }, finally = {
    plan(sequential)
  })

  inner_results[[j]] <- cell
}

profile <- list(
  index = i,
  scenario = scenario,
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
cat(sprintf("scenario=%d run=%d\n", scenario, run))
cat(sprintf("data generation + mice(m=50) imputation: %.1fs\n", profile$data_gen_elapsed))
cat("\ninner cells (model-fitting + bootstrap phase only; NA rows are cells where syrup itself failed):\n")
inner_df <- purrr::map_dfr(inner_results, function(x) {
  data.frame(workers1 = x$workers1, workers2 = x$workers2, grf_threads = x$grf_threads,
             CI_boot = x$CI_boot, ok = x$ok,
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
  cat(sprintf("\ntotal replicate walltime at each config (at this cell's CI_boot, not the real 200) = data_gen_elapsed + elapsed, e.g. best case %.1fs, worst case %.1fs\n",
              profile$data_gen_elapsed + min(inner_df$elapsed, na.rm = TRUE),
              profile$data_gen_elapsed + max(inner_df$elapsed, na.rm = TRUE)))
}

if (length(profile$warnings) > 0) {
  cat(sprintf("\n%d distinct warnings, most frequent:\n", length(profile$warnings)))
  print(head(sort(table(profile$warnings), decreasing = TRUE), 5))
}

output_dir <- file.path(dirname(path), "results", "missing", "ci_example", "profiling")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(profile, file.path(output_dir, paste0("prof_", i, ".RDS")))

print("Profiling run completed!")
