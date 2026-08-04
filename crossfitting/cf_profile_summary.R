##########
# title: turn the profiling runs into PBS directives for cf_1.sh
##########
# Reads the prof_*.RDS files written by cf_profile.R and the usage_*.txt logs
# written by the background sampler in cf_profile.sh, and prints the recommended
# walltime / mem / (workers, grf_threads) plus a time-by-arm table.

# libraries
library(dplyr)
library(tidyr)
library(purrr)
library(here)

# paths
path <- here()
prof_path <- file.path(dirname(path), "results", "crossfitting", "profiling")
log_path <- here("crossfitting", "jobscripts", "logs_profile")
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
         r_peak_mb = p$r_peak_mb, n_arms = p$n_arms)
})

cat("\n=== timings by configuration ===\n")
config <- runs %>%
  group_by(workers, grf_threads) %>%
  summarise(n_runs = n(),
            mean_elapsed = mean(elapsed),
            max_elapsed = max(elapsed),
            mean_cpu_seconds = mean(cpu_seconds),
            .groups = "drop") %>%
  arrange(mean_cpu_seconds)
print(config)

# ---- parallel efficiency ----------------------------------------------------
# t(1 worker) / (w * t(w workers)), at each thread count
cat("\n=== parallel efficiency ===\n")
efficiency <- config %>%
  group_by(grf_threads) %>%
  mutate(baseline = mean_elapsed[workers == 1],
         speedup = baseline / mean_elapsed,
         efficiency = speedup / workers) %>%
  ungroup() %>%
  select(workers, grf_threads, mean_elapsed, speedup, efficiency) %>%
  arrange(grf_threads, workers)
print(efficiency)

# the array runs throttled at %190, so it is throughput bound: pick the
# configuration that minimises cpu-seconds charged, not raw elapsed time
best <- config %>% slice_min(mean_cpu_seconds, n = 1)
cat(sprintf("\nchosen configuration: workers=%d, grf_threads=%d (fewest cpu-seconds)\n",
            best$workers, best$grf_threads))

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
print(arm_times, n = Inf)

# note: nuisance time is shared between the arms that use the same nuisance
# object, so summing mean_total over arms double counts. the replicate wall time
# in `config` above is the real per-job cost.

# ---- peak memory from the PBS accounting logs ------------------------------
# the background sampler in cf_profile.sh writes lines like
#   resources_used.mem = 3456789kb
parse_usage <- function(f) {
  txt <- readLines(f, warn = FALSE)
  mem <- regmatches(txt, regexpr("resources_used\\.mem = [0-9]+kb", txt))
  mem <- as.numeric(gsub("[^0-9]", "", mem))
  if (length(mem) == 0) return(NA_real_)
  max(mem) / 1024^2   # kb -> gb
}

usage_files <- list.files(log_path, pattern = "^usage_\\d+\\.txt$", full.names = TRUE)

if (length(usage_files) > 0) {
  peak_gb <- max(map_dbl(usage_files, parse_usage), na.rm = TRUE)
  mem_source <- "PBS accounting (includes multisession workers)"
} else {
  # fall back to the parent process's gc(), which misses the forked workers -
  # scale by the worker count as a crude stand-in and say so
  peak_gb <- max(runs$r_peak_mb) / 1024 * best$workers
  mem_source <- "R gc() x workers - FALLBACK, no usage_*.txt logs found"
  warning("no usage_*.txt logs found; the memory recommendation is a rough lower bound")
}

# ---- the recommendation -----------------------------------------------------
max_elapsed <- max(runs$elapsed[runs$workers == best$workers &
                                  runs$grf_threads == best$grf_threads])

walltime_sec <- ceiling(max_elapsed * walltime_factor / 1800) * 1800  # next 30 min
walltime_str <- sprintf("%02d:%02d:00", walltime_sec %/% 3600, (walltime_sec %% 3600) %/% 60)
mem_gb <- ceiling(peak_gb * mem_factor)
ncpus <- best$workers * best$grf_threads

cat("\n=== recommendation ===\n")
cat(sprintf("slowest profiled replicate : %.1f s\n", max_elapsed))
cat(sprintf("peak memory (%s): %.2f gb\n", mem_source, peak_gb))
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
