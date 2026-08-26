##########
# title: resource profiling for the interim-analysis validation array job
##########
# Runs one complete replicate - BOTH trial chunks, the way cts_val_analysis.R
# does - under a sweep of (interim_prop, workers, grf_threads), and records
# timings, memory and CPU so the PBS directives on cts_val_1.sh can be measured
# rather than guessed. Sister to continuous/cts_profile.R; see that file's
# header for the syrup/multisession rationale, not repeated here.
#
# ---- how to run it ----------------------------------------------------------
# Unlike continuous/cts_profile.R and binary/bin_profile.R, this is NOT a PBS
# array job - there is no cts_val_profile.sh. It is a loop, meant to be sourced
# inside an interactive RStudio session on the cluster:
#
#   request an RStudio session with 8 cores and 64gb, then
#   source(here::here("validation", "continuous", "cts_val_profile.R"))
#
# It writes one prof_<i>.RDS per cell as that cell finishes and skips cells
# already on disk, so an interrupted session is resumed by sourcing it again. A
# cell that errors is logged and skipped rather than ending the sweep. Set
# `cells` below to profile a subset - `cells <- 1` first is a good idea, to check
# syrup works in the session and to get a per-cell time before committing to the
# other 27.
#
# It has to be the cluster, not a laptop: syrup measures nothing on Windows, and
# a cell whose measurements are empty stops rather than saving a profile with
# -Inf memory in it, which would quietly poison the mem= recommendation.
#
# ---- the specific thing this settles ----------------------------------------
# cts_val_1.sh has always requested ncpus=5:ompthreads=5 to match
# cts_val_analysis.R's hardcoded `workers <- 5`, but num.threads was never
# passed to any grf call in cts_val_models.R - so every one of those 5
# multisession workers spawned forests on all visible cores regardless of what
# PBS allocated. cts_val_models.R and cts_val_analysis.R now accept/forward
# num.threads; this sweep measures whether the 5x1 shape was oversubscribed once
# that knob does something, and how few cores the job actually needs.
#
# There is a structural reason to expect `workers = 5` is the wrong shape. The
# only future_map() calls in this study are over the p covariates (the TE-VIM
# refits) and over the 9-combination xgboost tuning grid. nuisance_rf,
# run_causal_forest and run_dr_random_forest are a sequential spine that workers
# do nothing for and only num.threads speeds up. So the question is not "how
# many workers" but how to split a small core budget between the two.
#
# ---- what is and is not measured --------------------------------------------
# Both chunks are profiled, because both are what one array job pays for, and
# because the larger chunk drives peak memory. The rpart/lm block downstream of
# them in cts_val_analysis.R is not: it is sub-second against multi-minute
# forest fits, and the walltime factor in cts_val_profile_summary.R covers it.

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
source(here("R", "utils.R"))           # setup_rng_stream, timed
source(here("validation", "continuous", "cts_val_dgms.R"))
source(here("validation", "continuous", "cts_val_models.R"))
source(here("validation", "continuous", "cts_val_config.R"))

# ---- knobs ------------------------------------------------------------------
sample_interval <- 1        # seconds between syrup snapshots
overwrite <- FALSE          # TRUE re-runs cells that already have a prof_*.RDS

output_dir <- file.path(study$res_path, "profiling")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- the sweep --------------------------------------------------------------
# scenario and n are fixed at 3 and 1000 - unlike continuous/, this study's grid
# has only one of each (cts_val_config.R) - so interim_prop takes the place n
# holds in cts_profile.R's sweep. 0.75 is left out because it is the mirror of
# 0.25: the same pair of chunk sizes (250, 750) in the other order, so the same
# total cost.
available_cores <- as.integer(parallelly::availableCores())
if (available_cores < 2) {
  stop("this session sees ", available_cores, " core(s) - nothing to sweep. ",
       "request an RStudio session with 8.", call. = FALSE)
}
if (available_cores < 8) {
  warning("this session sees ", available_cores, " cores, not the 8 the sweep ",
          "is designed for - the larger configurations will be dropped and the ",
          "measured timings will reflect this session, not the job.",
          call. = FALSE, immediate. = TRUE)
}
max_cores <- min(8L, available_cores)

prof_params <- expand.grid(
  interim_prop = c(0.25, 0.5),
  run          = c(1, 2),
  workers      = c(1, 2, 4, 8),
  grf_threads  = c(1, 2),
  stringsAsFactors = FALSE
)
# a cell asking for more cores than the session has would measure contention,
# not the configuration - the same check cts_val_profile_summary.R makes again
# on the saved runs
prof_params <- prof_params[prof_params$workers * prof_params$grf_threads <= max_cores, ]
rownames(prof_params) <- NULL

# EDIT ME to profile a subset - e.g. cells <- 1 for a single trial cell.
# Trailing arguments override it, so the same file also runs as a quick
# non-interactive check: Rscript cts_val_profile.R 1
cells <- seq_len(nrow(prof_params))
cli_cells <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)))
if (length(cli_cells) > 0 && !anyNA(cli_cells)) cells <- cli_cells

cat(sprintf("sweep: %d cells over %d configurations, %d core(s) available\n",
            nrow(prof_params),
            nrow(distinct(prof_params, workers, grf_threads)),
            available_cores))

# ---- one cell ---------------------------------------------------------------
profile_cell <- function(i) {
  param <- prof_params[i, ]

  scenario <- 3
  n <- 1000
  interim_prop <- param$interim_prop
  run <- param$run
  workers <- param$workers
  grf_threads <- param$grf_threads

  # same seeding and n-dependent fold rule as cts_val_analysis.R, so profiled
  # cost matches production
  setup_rng_stream(run)

  gen1 <- generate_continuous_scenario_data(scenario, n * interim_prop)
  gen2 <- generate_continuous_scenario_data(scenario, n * (1 - interim_prop))

  data1 <- gen1$dataset
  data2 <- gen2$dataset

  n_folds1 <- ifelse(n * interim_prop < 250, 5, 10)
  n_folds2 <- ifelse(n * (1 - interim_prop) < 250, 5, 10)

  # the length of the TE-VIM future_map(): workers beyond this buy nothing on
  # that step, and the only other parallel step is a 9-combination grid
  p <- ncol(data1) - 2L

  Sys.setenv(OMP_NUM_THREADS = grf_threads)

  parent_pid <- Sys.getpid()

  # One R session serves every cell here, where the array-based profilers get a
  # fresh process per cell. The parent's RSS therefore ratchets upward across
  # cells, so the raw peak below overstates a real job's footprint by whatever
  # this session was already holding. Record the baseline and report both.
  # ps is syrup's own dependency, so it is here whenever syrup is, but a missing
  # baseline should degrade to NA rather than lose the cell's measurements
  baseline_rss_gb <- tryCatch(
    ps::ps_memory_info(ps::ps_handle())[["rss"]] / 1024^3,
    error = function(e) NA_real_
  )

  plan(multisession, workers = workers)

  # syrup() returns the measurement tibble, not the value of the expression - but
  # it evaluates the expression in this environment, so `total`, `res1` and
  # `res2` persist. verbose_timing = TRUE asks run_all_cate_methods() for a
  # per-block timing breakdown (results$timings), used below to build arm_times.
  usage <- syrup::syrup({
    total <- system.time({
      res1 <- run_all_cate_methods(data = data1, n_folds = n_folds1,
                                   num.threads = grf_threads,
                                   verbose_timing = TRUE)
      res2 <- run_all_cate_methods(data = data2, n_folds = n_folds2,
                                   num.threads = grf_threads,
                                   verbose_timing = TRUE)
    })
  }, interval = sample_interval)

  # explicit, not on.exit(): this is a loop at top level in an interactive
  # session, so on.exit would never fire and the previous cell's workers would
  # still be alive when the next one measures its process tree
  plan(sequential)

  # syrup reports every R process `ps` can see - in an RStudio session that
  # includes the IDE's own helpers - keep only this tree, as cts_profile.R does
  usage <- usage %>%
    filter(pid == parent_pid | ppid == parent_pid) %>%
    mutate(role = if_else(pid == parent_pid, "parent", "worker"),
           rss = as.numeric(rss),
           vms = as.numeric(vms))

  # An empty tibble here means syrup measured nothing at all - it is unreliable
  # on Windows, so a local smoke test lands in this branch. Stop rather than
  # carry on: every memory figure below would be -Inf, and saving that would put
  # a cell in the sweep that silently poisons the mem= recommendation.
  if (nrow(usage) == 0) {
    stop("cell ", i, ": syrup returned no measurements for this process tree, ",
         "so there is nothing to profile (the work itself took ",
         sprintf("%.1f", unname(total["elapsed"])), "s). syrup does not work on ",
         "Windows - run this on the cluster.", call. = FALSE)
  }

  # future falls back to sequential when workers == 1
  n_procs <- n_distinct(usage$pid)
  expected_procs <- if (workers > 1) workers + 1L else 1L
  if (n_procs != expected_procs) {
    warning(sprintf("cell %d: expected %d processes in the tree but saw %d",
                    i, expected_procs, n_procs))
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

  # per-block timings, one row per block run_all_cate_methods() timed, per chunk.
  # The chunk column is the point: it separates the sequential spine
  # (nuisance_rf, causal_forest, dr_random_forest) from the two parallel steps
  # (*_te_vims, *_shap_vims), which is what decides whether workers earn their
  # cores.
  arm_times <- bind_rows(
    tibble::tibble(chunk = 1L, arm = names(res1$timings),
                   time_total = unlist(res1$timings, use.names = FALSE)),
    tibble::tibble(chunk = 2L, arm = names(res2$timings),
                   time_total = unlist(res2$timings, use.names = FALSE))
  )

  profile <- list(
    index = i,
    scenario = scenario,
    n = n,
    interim_prop = interim_prop,
    n1 = nrow(data1),
    n2 = nrow(data2),
    p = p,
    run = run,
    workers = workers,
    grf_threads = grf_threads,
    n_folds1 = n_folds1,
    n_folds2 = n_folds2,
    elapsed = unname(total["elapsed"]),
    user = unname(total["user.self"]),
    system = unname(total["sys.self"]),
    # cpu-seconds charged: what PBS bills, and what 1100 unthrottled jobs are
    # really competing for
    cpu_seconds = unname(total["elapsed"]) * workers * grf_threads,
    arm_times = arm_times,
    usage = usage,                 # full time series, for the timeline plot
    by_snapshot = by_snapshot,
    by_process = by_process,
    n_procs = n_procs,
    expected_procs = expected_procs,
    peak_rss_gb = peak_rss_gb,
    baseline_rss_gb = baseline_rss_gb,
    peak_above_baseline_gb = peak_rss_gb - baseline_rss_gb,
    median_pct_cpu = median_pct_cpu,
    allocated_pct_cpu = allocated_pct_cpu,
    sample_interval = sample_interval,
    available_cores = available_cores,
    omp_num_threads = Sys.getenv("OMP_NUM_THREADS"),
    interactive_session = interactive(),
    r_version = R.version.string,
    timestamp = Sys.time()
  )

  saveRDS(profile, file.path(output_dir, paste0("prof_", i, ".RDS")))
  profile
}

# ---- the loop ---------------------------------------------------------------
elapsed_so_far <- numeric(0)

for (i in cells) {
  out_file <- file.path(output_dir, paste0("prof_", i, ".RDS"))
  if (file.exists(out_file) && !overwrite) {
    cat(sprintf("[%3d/%d] skipping - %s already exists\n",
                i, nrow(prof_params), basename(out_file)))
    next
  }

  param <- prof_params[i, ]
  cat(sprintf("\n[%3d/%d] interim_prop=%.2f run=%d workers=%d grf_threads=%d\n",
              i, nrow(prof_params), param$interim_prop, param$run,
              param$workers, param$grf_threads))

  # one bad cell should not end a sweep that is 27 cells from done - log it and
  # carry on, and sourcing the script again will retry it
  profile <- tryCatch(profile_cell(i), error = function(e) {
    plan(sequential)
    message(sprintf("        FAILED: %s", conditionMessage(e)))
    NULL
  })
  if (is.null(profile)) next

  elapsed_so_far <- c(elapsed_so_far, profile$elapsed)

  cat(sprintf("        elapsed=%.1fs  cpu-seconds=%.1f  peak=%.2fgb (%.2fgb above baseline)\n",
              profile$elapsed, profile$cpu_seconds, profile$peak_rss_gb,
              profile$peak_above_baseline_gb))
  cat(sprintf("        cpu while busy: %.0f%% observed against %.0f%% allocated%s\n",
              profile$median_pct_cpu, profile$allocated_pct_cpu,
              if (!is.na(profile$median_pct_cpu) &&
                  profile$median_pct_cpu > profile$allocated_pct_cpu * 1.1) {
                "  <-- OVERSUBSCRIBED"
              } else ""))

  # crude but honest: the remaining cells are not all the same cost (workers = 1
  # is the slow end), so this is a mean-based estimate, not a promise
  remaining <- cells[cells > i]
  # paste0() recycles a zero-length vector to "", so this filter has to be
  # guarded or an empty `remaining` comes back length 1 and the ETA claims a
  # cell that isn't there
  if (length(remaining) > 0) {
    remaining <- remaining[!file.exists(file.path(output_dir,
                                                  paste0("prof_", remaining, ".RDS")))]
  }
  if (length(remaining) > 0) {
    cat(sprintf("        ~%.0f min left over %d cells, at the mean so far\n",
                mean(elapsed_so_far) * length(remaining) / 60, length(remaining)))
  }

  # this session is reused by every cell, unlike the array's fresh process
  rm(profile)
  gc()
}

cat(sprintf("\n%d of %d cells on disk in %s\n",
            length(list.files(output_dir, pattern = "^prof_\\d+\\.RDS$")),
            nrow(prof_params), output_dir))
cat("next: source(here::here('validation', 'continuous', 'cts_val_profile_summary.R'))\n")
