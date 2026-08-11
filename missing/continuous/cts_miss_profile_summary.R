##########
# title: turn the profiling runs into PBS directives for cts_miss_1.sh/cts_miss_2.sh
##########
# Reads the prof_*.RDS files written by cts_miss_profile.R and prints the
# recommended walltime / mem / (workers, grf_threads), plus the
# observed-against-allocated CPU table, a parallel-efficiency table, a
# cost-by-method breakdown, and a time-by-model table. Everything comes from
# syrup's in-R measurements plus plain system.time() for the serial
# data-handling step, so there is no scheduler output to parse. Adapted from
# missing/binary/bin_miss_profile_summary.R; see that file for the parts that
# are unchanged, and cts_miss_profile.R's header for why this study needs an
# outer(method/scenario/mechanism/run)/inner(workers/grf_threads) split in
# the first place.
#
# Differences from bin_miss_profile_summary.R:
#
# 1. Each prof_<i>.RDS holds ONE outer row (method/scenario/mechanism/run)
#    with a list of 6 inner (workers, grf_threads) cells - flattened below
#    into one `runs` row per (index, workers, grf_threads), 72 x 6 = 432 rows
#    if the full sweep ran.
#
# 2. `total_elapsed = data_gen_elapsed + elapsed` is added per row: the
#    serial mice/missForest/etc step happens before plan(multisession) is
#    even set up, so it is additive to whichever (workers, grf_threads)
#    config handles the model-fitting phase, regardless of that config. The
#    walltime recommendation at the bottom sizes off `total_elapsed`, not
#    `elapsed` alone.
#
# 3. A "cost by method" table reports data_gen_elapsed pooled over
#    scenario/mechanism/run - the diagnostic this study specifically needs,
#    since cts_miss_config.R's MISS_METHODS_STUDY spans cheap methods
#    (complete_cases, mean_imputation, ...) and two expensive ones
#    (missforest's iterative imputation, multiple_imputation's mice(m = 50)).
#
# 4. TWO jobscripts get the directive rewrite, not one: cts_miss_1.sh (indices
#    1-5000) and cts_miss_2.sh (5001-9900) split the 9900-row production grid
#    between them, unlike missing/binary's single bin_miss_1.sh (1-9900). Both
#    tasks run the identical per-replicate workload, so the same computed
#    walltime/mem/select line and the same chosen (workers, grf_threads) are
#    written into both - only their #PBS -J ranges differ, and those are left
#    untouched.

# libraries
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(paletteer)
library(here)

# paths
path <- here()
prof_path <- file.path(dirname(path), "results", "missing", "continuous", "profiling")
jobscripts <- c(here("missing", "continuous", "jobscripts", "cts_miss_1.sh"),
                here("missing", "continuous", "jobscripts", "cts_miss_2.sh"))

# safety factors
walltime_factor <- 1.75   # a PBS job killed at the limit saves nothing
mem_factor <- 1.5

# ---- read the profiling runs -----------------------------------------------
prof_files <- list.files(prof_path, pattern = "^prof_\\d+\\.RDS$", full.names = TRUE)
stopifnot(length(prof_files) > 0)

profs <- map(prof_files, readRDS)

# one row per (index, workers, grf_threads) - the outer fields are repeated
# across each outer row's 6 inner cells. Cells where cts_miss_profile.R's
# tryCatch caught a syrup failure (cell$ok == FALSE, see that script's header
# on repeated syrup::syrup() calls per session) are dropped here rather than
# producing NA rows - they carry no usable timing/memory data.
skipped_cells <- 0L
runs <- map_dfr(profs, function(p) {
  map_dfr(p$inner, function(cell) {
    if (!isTRUE(cell$ok)) {
      skipped_cells <<- skipped_cells + 1L
      return(NULL)
    }
    tibble(index = p$index, method = p$method, scenario = p$scenario,
           mechanism = p$mechanism, run = p$run,
           data_gen_elapsed = p$data_gen_elapsed,
           workers = cell$workers, grf_threads = cell$grf_threads,
           elapsed = cell$elapsed, cpu_seconds = cell$cpu_seconds,
           peak_rss_gb = cell$peak_rss_gb, n_procs = cell$n_procs,
           expected_procs = cell$expected_procs,
           median_pct_cpu = cell$median_pct_cpu,
           allocated_pct_cpu = cell$allocated_pct_cpu,
           available_cores = p$available_cores, n_models = cell$n_models)
  })
}) %>%
  mutate(total_elapsed = data_gen_elapsed + elapsed)

if (skipped_cells > 0) {
  warning(sprintf("%d inner cells had no usable measurement (syrup failed on the profiling run) - excluded", skipped_cells))
}

# a cell profiled on a machine with fewer cores than it asked for is measuring
# that machine's contention, not the configuration
undersized <- runs %>% filter(available_cores < workers * grf_threads)
if (nrow(undersized) > 0) {
  warning(sprintf("%d profiling cells ran on a machine with fewer cores than the ",
                  nrow(undersized)),
          "configuration requested - their timings are not usable. ",
          "indices: ", paste(unique(undersized$index), collapse = ", "))
  print(as.data.frame(undersized[, c("index", "workers", "grf_threads", "available_cores")]),
        row.names = FALSE)
}

# future runs sequentially at workers == 1, so those cells have no child process
bad_procs <- runs %>% filter(n_procs != expected_procs)
if (nrow(bad_procs) > 0) {
  warning("some cells did not track the expected number of processes - ",
          "the memory figures for those may be wrong. indices: ",
          paste(unique(bad_procs$index), collapse = ", "))
}

# ---- cost by method (the new axis this study introduces) -------------------
cat("\n=== data generation + missingness handling cost by method (pooled over scenario/mechanism/run) ===\n")
# one row per outer index first, so 6 repeated inner cells don't inflate the mean
by_method <- runs %>%
  distinct(index, method, data_gen_elapsed) %>%
  group_by(method) %>%
  summarise(n_runs = n(),
            mean_data_gen_s = mean(data_gen_elapsed),
            max_data_gen_s = max(data_gen_elapsed),
            .groups = "drop") %>%
  arrange(desc(mean_data_gen_s))
print(as.data.frame(by_method), digits = 3, row.names = FALSE)

cat("\n=== timings by configuration (pooled over method/scenario/mechanism/run; model-fitting phase only) ===\n")
config <- runs %>%
  group_by(workers, grf_threads) %>%
  summarise(n_runs = n(),
            mean_elapsed = mean(elapsed),
            max_elapsed = max(elapsed),
            mean_cpu_seconds = mean(cpu_seconds),
            mean_peak_rss_gb = mean(peak_rss_gb),
            max_peak_rss_gb = max(peak_rss_gb),
            .groups = "drop") %>%
  arrange(mean_cpu_seconds)
print(as.data.frame(config), digits = 3, row.names = FALSE)

# ---- oversubscription -------------------------------------------------------
# the question the sweep exists to answer, measured rather than inferred
cat("\n=== CPU: observed against allocated ===\n")
cpu <- runs %>%
  group_by(workers, grf_threads) %>%
  summarise(allocated_pct = first(allocated_pct_cpu),
            observed_pct = median(median_pct_cpu, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(ratio = observed_pct / allocated_pct,
         verdict = case_when(ratio > 1.1 ~ "oversubscribed",
                             ratio < 0.6 ~ "under-used",
                             TRUE ~ "ok")) %>%
  arrange(workers, grf_threads)
print(as.data.frame(cpu), digits = 3, row.names = FALSE)

# ---- parallel efficiency ----------------------------------------------------
cat("\n=== parallel efficiency ===\n")
# the workers == 1 cell is the baseline. if it failed for a thread count, that
# whole row would subset to numeric(0) and error, so fall back to NA instead
efficiency <- config %>%
  group_by(grf_threads) %>%
  mutate(baseline = if (any(workers == 1)) mean_elapsed[workers == 1][1] else NA_real_,
         speedup = baseline / mean_elapsed,
         efficiency = speedup / workers) %>%
  ungroup() %>%
  select(workers, grf_threads, mean_elapsed, speedup, efficiency) %>%
  arrange(grf_threads, workers)
if (anyNA(efficiency$speedup)) {
  warning("no workers == 1 cell for at least one thread count - ",
          "efficiency is NA there")
}
print(as.data.frame(efficiency), digits = 3, row.names = FALSE)

# the array runs throttled at %190 (same as binary's bin_miss_1.sh), so it is
# throughput bound: pick the configuration that minimises mean cpu-seconds
# charged across the whole sweep (all methods/scenarios/mechanisms pooled,
# since one PBS resource line serves every arm in both array jobs)
best <- config %>% slice_min(mean_cpu_seconds, n = 1)
cat(sprintf("\nchosen configuration: workers=%d, grf_threads=%d (fewest cpu-seconds)\n",
            best$workers, best$grf_threads))

best_runs <- runs %>% filter(workers == best$workers, grf_threads == best$grf_threads)

# ---- where the time actually goes ------------------------------------------
cat("\n=== mean time by model, at the chosen configuration (model-fitting phase only) ===\n")
arm_times <- map_dfr(profs, function(p) {
  cell <- keep(p$inner, ~ isTRUE(.x$ok) && .x$workers == best$workers && .x$grf_threads == best$grf_threads)
  if (length(cell) == 0) return(NULL)
  mutate(cell[[1]]$arm_times, method = p$method)
}) %>%
  group_by(arm) %>%
  summarise(mean_time = mean(time_total),
            max_time = max(time_total),
            .groups = "drop") %>%
  arrange(desc(mean_time))
print(as.data.frame(arm_times), digits = 3, row.names = FALSE)

# note: nuisance_rf is shared overhead computed once and reused by several
# models below it, not a per-model cost in its own right - summing mean_time
# over rows double counts it. multiple_imputation's arm_times are a 50-way
# SUM (cts_miss_profile.R's header explains why), so they are not comparable
# 1:1 against the other methods' single-call times - treat them as "total
# model-fitting CPU time this method needs per replicate", not "time per model
# fit". The replicate wall time in `config`/the recommendation below is the
# real per-job cost either way.

# ---- per process behaviour --------------------------------------------------
cat("\n=== per process, at the chosen configuration ===\n")
by_process <- map_dfr(profs, function(p) {
  cell <- keep(p$inner, ~ isTRUE(.x$ok) && .x$workers == best$workers && .x$grf_threads == best$grf_threads)
  if (length(cell) == 0) return(NULL)
  mutate(cell[[1]]$by_process, index = p$index)
}) %>%
  group_by(role) %>%
  summarise(n = n(),
            mean_max_rss_gb = mean(max_rss_gb),
            mean_median_pct_cpu = mean(median_pct_cpu),
            mean_max_pct_cpu = mean(max_pct_cpu),
            .groups = "drop")
print(as.data.frame(by_process), digits = 3, row.names = FALSE)

# ---- timeline plot for one representative cell ------------------------------
# the slowest replicate (data-gen + model-fitting) at the chosen
# configuration - almost always a multiple_imputation or missforest one,
# which is exactly the case that should dominate the sizing below.
# Wrapped in tryCatch: the recommendation and jobscript rewrite below are the
# actual deliverable of this script, so a plotting failure (e.g. a cell whose
# usage tibble came back empty - see cts_miss_profile.R's header on syrup
# occasionally failing mid-sweep) must not stop the script before it gets there.
tryCatch({
  rep_index <- best_runs$index[which.max(best_runs$total_elapsed)]
  rep_prof <- profs[[which(map_dbl(profs, "index") == rep_index)]]
  rep_cell <- keep(rep_prof$inner, ~ isTRUE(.x$ok) && .x$workers == best$workers && .x$grf_threads == best$grf_threads)[[1]]

  timeline <- rep_cell$usage %>%
    mutate(secs = as.numeric(difftime(time, min(time), units = "secs")),
           rss_gb = rss / 1024^3) %>%
    select(secs, pid, role, rss_gb, pct_cpu) %>%
    pivot_longer(c(rss_gb, pct_cpu), names_to = "measure", values_to = "value") %>%
    mutate(measure = recode(measure, rss_gb = "Memory (gb)", pct_cpu = "CPU (%)"))

  timeline_plot <- timeline %>%
    ggplot(aes(x = secs, y = value, colour = role, group = factor(pid))) +
    geom_line(alpha = 0.8) +
    facet_wrap(~measure, ncol = 1, scales = "free_y") +
    scale_colour_paletteer_d("rcartocolor::Safe") +
    theme_minimal() +
    labs(title = sprintf("Resource timeline: workers=%d, grf_threads=%d, method=%s, scenario %d, mechanism %s",
                         best$workers, best$grf_threads, rep_prof$method,
                         rep_prof$scenario, rep_prof$mechanism),
         subtitle = "One line per process; sampled by syrup. Data generation/handling (not shown) happens before this window.",
         x = "Seconds since model-fitting started", y = NULL, colour = NULL)
  ggsave("cts_miss_profile_timeline.png", plot = timeline_plot, path = prof_path,
         width = 21, height = 15, units = "cm")
  cat(sprintf("\ntimeline plot written to %s\n", file.path(prof_path, "cts_miss_profile_timeline.png")))
}, error = function(e) {
  warning("timeline plot failed, skipping it - the recommendation below is unaffected: ",
          conditionMessage(e))
})

# ---- the recommendation -----------------------------------------------------
# total_elapsed (data-gen + model-fitting), not elapsed alone, since the
# serial data-handling step happens inside the same PBS walltime budget
max_elapsed <- max(best_runs$total_elapsed)
peak_gb <- max(best_runs$peak_rss_gb)
ncpus <- best$workers * best$grf_threads

# peak_gb can come back -Inf if every cell at the chosen configuration had an
# empty usage tibble (the local-Windows syrup process-tracking gap
# cts_miss_profile.R's header describes - on a real HPC node with a working ps,
# this shouldn't happen). Refuse to size or write a nonsensical mem= from that
# rather than let sprintf("%d", -Inf) error out below the walltime is already
# printed, or silently write "mem=-Infgb" into the real jobscripts.
if (!is.finite(peak_gb)) {
  warning("peak_rss_gb was not finite for every cell at the chosen configuration - ",
          "no usable memory measurement (see the process-count warnings above). ",
          "Skipping the mem= recommendation and the jobscript rewrite; rerun the ",
          "profiling sweep somewhere syrup can see the process tree (e.g. the HPC login node).")
  mem_gb <- NA_real_
} else {
  mem_gb <- ceiling(peak_gb * mem_factor)
}

walltime_sec <- ceiling(max_elapsed * walltime_factor / 1800) * 1800  # next 30 min
walltime_str <- sprintf("%02d:%02d:00", walltime_sec %/% 3600, (walltime_sec %% 3600) %/% 60)

cat("\n=== recommendation ===\n")
cat(sprintf("slowest profiled replicate (data-gen + model-fitting): %.1f s\n", max_elapsed))
cat(sprintf("peak memory (summed RSS across the process tree, model-fitting phase): %.2f gb\n", peak_gb))
cat("  this overcounts shared library pages, so it is an upper bound - the safe\n")
cat("  direction for a mem= request. cross-check once against qstat -fx on the\n")
cat("  first real subjobs; PBS's cgroup figure should sit below it.\n")
cat(sprintf("\n#PBS -l walltime=%s\n", walltime_str))
if (is.na(mem_gb)) {
  cat(sprintf("#PBS -l select=1:ncpus=%d:ompthreads=%d:mem=<no usable measurement>\n", ncpus, ncpus))
} else {
  cat(sprintf("#PBS -l select=1:ncpus=%d:ompthreads=%d:mem=%dgb\n", ncpus, ncpus, mem_gb))
}
cat(sprintf("Rscript cts_miss_analysis.R \"$PBS_ARRAY_INDEX\" %d %d\n\n",
            best$workers, best$grf_threads))

# ---- write the directives into cts_miss_1.sh and cts_miss_2.sh --------------
# workers/grf_threads are written as trailing args on the Rscript line, not as
# a manual edit to cts_miss_analysis.R's defaults, so the PBS resource request
# and the R-level parallelism config can never drift apart. Both jobscripts
# get the identical rewrite (see this file's header, point 4) - only their
# #PBS -J ranges differ, and those are left untouched. Skipped entirely if
# mem_gb couldn't be sized - a walltime-only rewrite would leave a stale mem=
# untouched next to a fresh walltime, silently mismatched.
if (is.na(mem_gb)) {
  cat("jobscripts NOT rewritten - no usable memory measurement (see warning above)\n")
} else {
  for (jobscript in jobscripts) {
    if (file.exists(jobscript)) {
      lines <- readLines(jobscript)
      lines <- sub("^#PBS -l walltime=.*$",
                   sprintf("#PBS -l walltime=%s", walltime_str), lines)
      lines <- sub("^#PBS -l select=.*$",
                   sprintf("#PBS -l select=1:ncpus=%d:ompthreads=%d:mem=%dgb",
                           ncpus, ncpus, mem_gb), lines)
      lines <- sub('^Rscript cts_miss_analysis\\.R "\\$PBS_ARRAY_INDEX".*$',
                   sprintf('Rscript cts_miss_analysis.R "$PBS_ARRAY_INDEX" %d %d',
                           best$workers, best$grf_threads),
                   lines)
      writeLines(lines, jobscript)
      cat(sprintf("directives and workers/grf_threads written into %s\n", jobscript))
    } else {
      warning(sprintf("%s not found - directives printed only", jobscript))
    }
  }
}
