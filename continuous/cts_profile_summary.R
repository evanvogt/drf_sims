##########
# title: turn the profiling runs into PBS directives for cts_1.sh
##########
# Reads the prof_*.RDS files written by cts_profile.R and prints the recommended
# walltime / mem / (workers, grf_threads), plus the observed-against-allocated CPU
# table, a parallel-efficiency table, a per-n cost breakdown, and a time-by-model
# table. Everything comes from syrup's in-R measurements, so there is no
# scheduler output to parse. Direct port of crossfitting/cf_profile_summary.R;
# see that file for the parts that are unchanged.
#
# One difference from the crossfitting version: continuous's single array job
# spans n in {100, 250, 500, 1000} under one #PBS -l line (crossfitting fixes
# n = 500), so `n` is carried through the profiling tibble and printed as its
# own breakdown - but the (workers, grf_threads) choice and the walltime/mem
# sizing below still pool across n, same as they pool across scenario/run,
# because one PBS resource line has to cover the whole array either way. The
# max()-based sizing at the bottom naturally lands on the n = 1000 cell since
# that is the slowest/heaviest one profiled.

# libraries
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(paletteer)
library(here)

# paths
path <- here()
prof_path <- file.path(dirname(path), "results", "continuous", "profiling")
jobscript <- here("continuous", "jobscripts", "cts_1.sh")

# safety factors
walltime_factor <- 1.75   # a PBS job killed at the limit saves nothing
mem_factor <- 1.5

# ---- read the profiling runs -----------------------------------------------
prof_files <- list.files(prof_path, pattern = "^prof_\\d+\\.RDS$", full.names = TRUE)
stopifnot(length(prof_files) > 0)

profs <- map(prof_files, readRDS)

runs <- map_dfr(profs, function(p) {
  tibble(index = p$index, scenario = p$scenario, n = p$n, run = p$run,
         workers = p$workers, grf_threads = p$grf_threads,
         elapsed = p$elapsed, cpu_seconds = p$cpu_seconds,
         peak_rss_gb = p$peak_rss_gb, n_procs = p$n_procs,
         median_pct_cpu = p$median_pct_cpu, allocated_pct_cpu = p$allocated_pct_cpu,
         available_cores = p$available_cores, n_models = p$n_models)
})

# a cell profiled on a machine with fewer cores than it asked for is measuring
# that machine's contention, not the configuration
undersized <- runs %>% filter(available_cores < workers * grf_threads)
if (nrow(undersized) > 0) {
  warning(sprintf("%d profiling cells ran on a machine with fewer cores than the ",
                  nrow(undersized)),
          "configuration requested - their timings are not usable. ",
          "indices: ", paste(undersized$index, collapse = ", "))
  print(as.data.frame(undersized[, c("index", "workers", "grf_threads", "available_cores")]),
        row.names = FALSE)
}

# future runs sequentially at workers == 1, so those cells have no child process
bad_procs <- runs %>% mutate(expected = if_else(workers > 1, workers + 1L, 1L)) %>%
  filter(n_procs != expected)
if (nrow(bad_procs) > 0) {
  warning("some cells did not track the expected number of processes - ",
          "the memory figures for those may be wrong. indices: ",
          paste(bad_procs$index, collapse = ", "))
}

# ---- cost by n ---------------------------------------------------------------
# continuous-specific: how much n alone drives cost, pooling over the other axes
cat("\n=== timings by n (pooled over scenario/workers/grf_threads/run) ===\n")
by_n <- runs %>%
  group_by(n) %>%
  summarise(n_runs = n(),
            mean_elapsed = mean(elapsed),
            max_elapsed = max(elapsed),
            mean_peak_rss_gb = mean(peak_rss_gb),
            max_peak_rss_gb = max(peak_rss_gb),
            .groups = "drop") %>%
  arrange(n)
print(as.data.frame(by_n), digits = 3, row.names = FALSE)

cat("\n=== timings by configuration (pooled over scenario/n/run) ===\n")
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

# the array runs throttled at %190, so it is throughput bound: pick the
# configuration that minimises mean cpu-seconds charged across the whole sweep
# (all n pooled together, since one PBS resource line serves every n in the array)
best <- config %>% slice_min(mean_cpu_seconds, n = 1)
cat(sprintf("\nchosen configuration: workers=%d, grf_threads=%d (fewest cpu-seconds)\n",
            best$workers, best$grf_threads))

best_runs <- runs %>% filter(workers == best$workers, grf_threads == best$grf_threads)

# ---- where the time actually goes ------------------------------------------
cat("\n=== mean time by model, at the chosen configuration ===\n")
arm_times <- map_dfr(profs, function(p) {
  if (p$workers != best$workers || p$grf_threads != best$grf_threads) return(NULL)
  p$arm_times
}) %>%
  group_by(arm) %>%
  summarise(mean_time = mean(time_total),
            max_time = max(time_total),
            .groups = "drop") %>%
  arrange(desc(mean_time))
print(as.data.frame(arm_times), digits = 3, row.names = FALSE)

# note: nuisance_rf is shared overhead computed once and reused by several
# models below it, not a per-model cost in its own right - summing mean_time
# over rows double counts it. The replicate wall time in `config` above is the
# real per-job cost.

# ---- per process behaviour --------------------------------------------------
cat("\n=== per process, at the chosen configuration ===\n")
by_process <- map_dfr(profs, function(p) {
  if (p$workers != best$workers || p$grf_threads != best$grf_threads) return(NULL)
  mutate(p$by_process, index = p$index)
}) %>%
  group_by(role) %>%
  summarise(n = n(),
            mean_max_rss_gb = mean(max_rss_gb),
            mean_median_pct_cpu = mean(median_pct_cpu),
            mean_max_pct_cpu = mean(max_pct_cpu),
            .groups = "drop")
print(as.data.frame(by_process), digits = 3, row.names = FALSE)

# ---- timeline plot for one representative cell ------------------------------
# the slowest cell at the chosen configuration - almost always an n = 1000 one,
# which is exactly the case that should dominate the sizing below
rep_prof <- profs[[which.max(vapply(profs, function(p) {
  if (p$workers != best$workers || p$grf_threads != best$grf_threads) return(-1)
  p$elapsed
}, numeric(1)))]]

timeline <- rep_prof$usage %>%
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
  labs(title = sprintf("Resource timeline: workers=%d, grf_threads=%d, scenario %d, n=%d",
                       rep_prof$workers, rep_prof$grf_threads, rep_prof$scenario, rep_prof$n),
       subtitle = "One line per process; sampled by syrup",
       x = "Seconds since start", y = NULL, colour = NULL)
ggsave("cts_profile_timeline.png", plot = timeline_plot, path = prof_path,
       width = 21, height = 15, units = "cm")
cat(sprintf("\ntimeline plot written to %s\n", file.path(prof_path, "cts_profile_timeline.png")))

# ---- the recommendation -----------------------------------------------------
max_elapsed <- max(best_runs$elapsed)
peak_gb <- max(best_runs$peak_rss_gb)

walltime_sec <- ceiling(max_elapsed * walltime_factor / 1800) * 1800  # next 30 min
walltime_str <- sprintf("%02d:%02d:00", walltime_sec %/% 3600, (walltime_sec %% 3600) %/% 60)
mem_gb <- ceiling(peak_gb * mem_factor)
ncpus <- best$workers * best$grf_threads

cat("\n=== recommendation ===\n")
cat(sprintf("slowest profiled replicate : %.1f s\n", max_elapsed))
cat(sprintf("peak memory (summed RSS across the process tree): %.2f gb\n", peak_gb))
cat("  this overcounts shared library pages, so it is an upper bound - the safe\n")
cat("  direction for a mem= request. cross-check once against qstat -fx on the\n")
cat("  first real subjobs; PBS's cgroup figure should sit below it.\n")
cat(sprintf("\n#PBS -l walltime=%s\n", walltime_str))
cat(sprintf("#PBS -l select=1:ncpus=%d:ompthreads=%d:mem=%dgb\n", ncpus, ncpus, mem_gb))
cat(sprintf("Rscript cts_analysis.R \"$PBS_ARRAY_INDEX\" %d %d\n\n",
            best$workers, best$grf_threads))

# ---- write the directives into cts_1.sh -------------------------------------
# workers/grf_threads are written as trailing args on the Rscript line, not as a
# manual edit to cts_analysis.R's defaults, so the PBS resource request and the
# R-level parallelism config can never drift apart.
if (file.exists(jobscript)) {
  lines <- readLines(jobscript)
  lines <- sub("^#PBS -l walltime=.*$",
               sprintf("#PBS -l walltime=%s", walltime_str), lines)
  lines <- sub("^#PBS -l select=.*$",
               sprintf("#PBS -l select=1:ncpus=%d:ompthreads=%d:mem=%dgb",
                       ncpus, ncpus, mem_gb), lines)
  lines <- sub('^Rscript cts_analysis\\.R "\\$PBS_ARRAY_INDEX".*$',
               sprintf('Rscript cts_analysis.R "$PBS_ARRAY_INDEX" %d %d',
                       best$workers, best$grf_threads),
               lines)
  writeLines(lines, jobscript)
  cat(sprintf("directives and workers/grf_threads written into %s\n", jobscript))
} else {
  warning("cts_1.sh not found - directives printed only")
}
