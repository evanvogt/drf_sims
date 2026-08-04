##########
# title: turn the profiling runs into PBS directives for cf_1.sh
##########
# Reads the prof_*.RDS files written by cf_profile.R and prints the recommended
# walltime / mem / (workers, grf_threads), plus the observed-against-allocated CPU
# table and a time-by-arm table. Everything comes from syrup's in-R measurements,
# so there is no scheduler output to parse.

# libraries
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(paletteer)
library(here)

# paths
path <- here()
prof_path <- file.path(dirname(path), "results", "crossfitting", "profiling")
jobscript <- here("crossfitting", "jobscripts", "cf_1.sh")

# safety factors
walltime_factor <- 1.75   # a PBS job killed at the limit saves nothing
mem_factor <- 1.5

# ---- read the profiling runs -----------------------------------------------
prof_files <- list.files(prof_path, pattern = "^prof_\\d+\\.RDS$", full.names = TRUE)
stopifnot(length(prof_files) > 0)

profs <- map(prof_files, readRDS)

runs <- map_dfr(profs, function(p) {
  tibble(index = p$index, scenario = p$scenario, run = p$run,
         workers = p$workers, grf_threads = p$grf_threads,
         elapsed = p$elapsed, cpu_seconds = p$cpu_seconds,
         peak_rss_gb = p$peak_rss_gb, n_procs = p$n_procs,
         median_pct_cpu = p$median_pct_cpu, allocated_pct_cpu = p$allocated_pct_cpu,
         available_cores = p$available_cores, n_arms = p$n_arms)
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

bad_procs <- runs %>% filter(n_procs != workers + 1)
if (nrow(bad_procs) > 0) {
  warning("some cells did not track the expected number of processes - ",
          "the memory figures for those may be wrong. indices: ",
          paste(bad_procs$index, collapse = ", "))
}

cat("\n=== timings by configuration ===\n")
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
efficiency <- config %>%
  group_by(grf_threads) %>%
  mutate(baseline = mean_elapsed[workers == 1],
         speedup = baseline / mean_elapsed,
         efficiency = speedup / workers) %>%
  ungroup() %>%
  select(workers, grf_threads, mean_elapsed, speedup, efficiency) %>%
  arrange(grf_threads, workers)
print(as.data.frame(efficiency), digits = 3, row.names = FALSE)

# the array runs throttled at %190, so it is throughput bound: pick the
# configuration that minimises cpu-seconds charged, not raw elapsed time
best <- config %>% slice_min(mean_cpu_seconds, n = 1)
cat(sprintf("\nchosen configuration: workers=%d, grf_threads=%d (fewest cpu-seconds)\n",
            best$workers, best$grf_threads))

best_runs <- runs %>% filter(workers == best$workers, grf_threads == best$grf_threads)

# ---- where the time actually goes ------------------------------------------
cat("\n=== mean time by arm, at the chosen configuration ===\n")
arm_times <- map_dfr(profs, function(p) {
  if (p$workers != best$workers || p$grf_threads != best$grf_threads) return(NULL)
  p$arm_times
}) %>%
  group_by(arm, family, variant) %>%
  summarise(mean_nuisance = mean(time_nuisance),
            mean_stage2 = mean(time_stage2),
            mean_total = mean(time_total),
            .groups = "drop") %>%
  arrange(desc(mean_total))
print(as.data.frame(arm_times), digits = 3, row.names = FALSE)

# note: nuisance time is shared between the arms that use the same nuisance
# object, so summing mean_total over arms double counts. the replicate wall time
# in `config` above is the real per-job cost.

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
# stage boundaries are visible in this, which is the quickest way to see whether
# one nuisance computation dominates memory or leaves the workers idle
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
  labs(title = sprintf("Resource timeline: workers=%d, grf_threads=%d, scenario %d",
                       rep_prof$workers, rep_prof$grf_threads, rep_prof$scenario),
       subtitle = "One line per process; sampled by syrup",
       x = "Seconds since start", y = NULL, colour = NULL)
ggsave("cf_profile_timeline.png", plot = timeline_plot, path = prof_path,
       width = 21, height = 15, units = "cm")
cat(sprintf("\ntimeline plot written to %s\n", file.path(prof_path, "cf_profile_timeline.png")))

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
cat(sprintf("#PBS -l select=1:ncpus=%d:ompthreads=%d:mem=%dgb\n", ncpus, best$grf_threads, mem_gb))
cat(sprintf("\nand in cf_analysis.R: workers <- %d ; grf_threads <- %d\n\n",
            best$workers, best$grf_threads))

# ---- write the directives into cf_1.sh -------------------------------------
if (file.exists(jobscript)) {
  lines <- readLines(jobscript)
  lines <- sub("^#PBS -l walltime=.*$",
               sprintf("#PBS -l walltime=%s", walltime_str), lines)
  lines <- sub("^#PBS -l select=.*$",
               sprintf("#PBS -l select=1:ncpus=%d:ompthreads=%d:mem=%dgb",
                       ncpus, best$grf_threads, mem_gb), lines)
  writeLines(lines, jobscript)
  cat(sprintf("directives written into %s\n", jobscript))
  cat("check that workers / grf_threads in cf_analysis.R match before submitting.\n")
} else {
  warning("cf_1.sh not found - directives printed only")
}
