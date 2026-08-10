##########
# title: turn the profiling runs into PBS directives for me_1.sh
##########
# Reads the prof_*.RDS files written by me_profile.R and prints the
# recommended walltime / mem / (workers, n_cores), plus the
# observed-against-allocated CPU table and a phase-timing table (candidate
# models vs. nuisance-evaluation pipelines - this study's analogue of
# cf_profile_summary.R's per-arm timing table, since candidate models here
# don't carry individual timings the way crossfitting/'s arms do).
# Everything comes from syrup's in-R measurements, so there is no scheduler
# output to parse.

# libraries
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(paletteer)
library(here)

# paths
path <- here()
prof_path <- file.path(dirname(path), "results", "model_evaluation", "profiling")
jobscript <- here("model_evaluation", "jobscripts", "me_1.sh")

# safety factors - larger than crossfitting/cf_profile_summary.R's (1.75,
# 1.5): H2O AutoML has no max_runtime_secs cap in this study (see README.md),
# so per-replicate cost is inherently less predictable than a fixed set of
# RF/causal-forest fits, and a killed PBS job saves nothing.
walltime_factor <- 2.5
mem_factor <- 1.5

# ---- read the profiling runs -----------------------------------------------
prof_files <- list.files(prof_path, pattern = "^prof_\\d+\\.RDS$", full.names = TRUE)
stopifnot(length(prof_files) > 0)

profs <- map(prof_files, readRDS)

runs <- map_dfr(profs, function(p) {
  tibble(index = p$index, scenario = p$scenario, run = p$run,
         workers = p$workers, n_cores = p$n_cores,
         elapsed_models = p$elapsed_models, elapsed_nuisances = p$elapsed_nuisances,
         elapsed = p$elapsed, cpu_seconds = p$cpu_seconds,
         peak_rss_gb = p$peak_rss_gb, n_procs = p$n_procs,
         median_pct_cpu = p$median_pct_cpu, allocated_pct_cpu = p$allocated_pct_cpu,
         available_cores = p$available_cores)
})

# a cell profiled on a machine with fewer cores than it asked for is
# measuring that machine's contention, not the configuration. workers and
# n_cores parallelise sequential, not concurrent, phases - so the relevant
# ceiling is max(workers, n_cores), not their sum or product.
undersized <- runs %>% filter(available_cores < pmax(workers, n_cores))
if (nrow(undersized) > 0) {
  warning(sprintf("%d profiling cells ran on a machine with fewer cores than the ",
                  nrow(undersized)),
          "configuration requested - their timings are not usable. ",
          "indices: ", paste(undersized$index, collapse = ", "))
  print(as.data.frame(undersized[, c("index", "workers", "n_cores", "available_cores")]),
        row.names = FALSE)
}

bad_procs <- runs %>% mutate(expected = workers + 1L) %>%
  filter(n_procs != expected)
if (nrow(bad_procs) > 0) {
  warning("some cells did not track the expected number of R processes - the ",
          "memory figures for those may be wrong (and may be missing H2O's ",
          "JVM entirely - see me_profile.R's header). indices: ",
          paste(bad_procs$index, collapse = ", "))
}

cat("\n=== timings by configuration ===\n")
config <- runs %>%
  group_by(workers, n_cores) %>%
  summarise(n_runs = n(),
            mean_elapsed_models = mean(elapsed_models),
            mean_elapsed_nuisances = mean(elapsed_nuisances),
            mean_elapsed = mean(elapsed),
            max_elapsed = max(elapsed),
            mean_cpu_seconds = mean(cpu_seconds),
            mean_peak_rss_gb = mean(peak_rss_gb),
            max_peak_rss_gb = max(peak_rss_gb),
            .groups = "drop") %>%
  arrange(mean_cpu_seconds)
print(as.data.frame(config), digits = 3, row.names = FALSE)

# ---- where the time actually goes -------------------------------------------
cat("\n=== phase breakdown: candidate models vs. nuisance-evaluation pipelines ===\n")
cat("(this is model_evaluation's analogue of cf_profile_summary.R's per-arm\n")
cat(" table - individual candidate models aren't separately timed here)\n")
phase <- config %>%
  mutate(pct_models = 100 * mean_elapsed_models / mean_elapsed,
         pct_nuisances = 100 * mean_elapsed_nuisances / mean_elapsed) %>%
  select(workers, n_cores, mean_elapsed_models, mean_elapsed_nuisances, pct_nuisances)
print(as.data.frame(phase), digits = 3, row.names = FALSE)

# ---- oversubscription -------------------------------------------------------
cat("\n=== CPU: observed against allocated (approximate, see me_profile.R note) ===\n")
cpu <- runs %>%
  group_by(workers, n_cores) %>%
  summarise(allocated_pct = first(allocated_pct_cpu),
            observed_pct = median(median_pct_cpu, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(ratio = observed_pct / allocated_pct,
         verdict = case_when(ratio > 1.1 ~ "oversubscribed",
                             ratio < 0.6 ~ "under-used",
                             TRUE ~ "ok")) %>%
  arrange(workers, n_cores)
print(as.data.frame(cpu), digits = 3, row.names = FALSE)

# the array will be throttled (H2O's per-task memory footprint rules out
# anything like continuous's %190 or crossfitting's %380 - see README.md), so
# it is throughput bound: pick the configuration that minimises cpu-seconds
# charged, not raw elapsed time.
best <- config %>% slice_min(mean_cpu_seconds, n = 1)
cat(sprintf("\nchosen configuration: workers=%d, n_cores=%d (fewest cpu-seconds)\n",
            best$workers, best$n_cores))

best_runs <- runs %>% filter(workers == best$workers, n_cores == best$n_cores)

# ---- per process behaviour --------------------------------------------------
cat("\n=== per process, at the chosen configuration ===\n")
by_process <- map_dfr(profs, function(p) {
  if (p$workers != best$workers || p$n_cores != best$n_cores) return(NULL)
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
rep_prof <- profs[[which.max(vapply(profs, function(p) {
  if (p$workers != best$workers || p$n_cores != best$n_cores) return(-1)
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
  geom_vline(xintercept = rep_prof$elapsed_models, linetype = "dashed", alpha = 0.4) +
  facet_wrap(~measure, ncol = 1, scales = "free_y") +
  scale_colour_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(title = sprintf("Resource timeline: workers=%d, n_cores=%d, scenario %d",
                       rep_prof$workers, rep_prof$n_cores, rep_prof$scenario),
       subtitle = "One line per process; dashed line marks the models->nuisances phase boundary",
       x = "Seconds since start", y = NULL, colour = NULL)
ggsave("me_profile_timeline.png", plot = timeline_plot, path = prof_path,
       width = 21, height = 15, units = "cm")
cat(sprintf("\ntimeline plot written to %s\n", file.path(prof_path, "me_profile_timeline.png")))

# ---- the recommendation -----------------------------------------------------
max_elapsed <- max(best_runs$elapsed)
peak_gb <- max(best_runs$peak_rss_gb)

walltime_sec <- ceiling(max_elapsed * walltime_factor / 1800) * 1800 # next 30 min
walltime_str <- sprintf("%02d:%02d:00", walltime_sec %/% 3600, (walltime_sec %% 3600) %/% 60)
mem_gb <- ceiling(peak_gb * mem_factor)
ncpus <- max(best$workers, best$n_cores) # sequential phases, not additive

cat("\n=== recommendation ===\n")
cat(sprintf("slowest profiled replicate : %.1f s\n", max_elapsed))
cat(sprintf("peak memory (summed RSS across the tracked R process tree): %.2f gb\n", peak_gb))
cat("  this overcounts shared library pages (safe direction for mem=), UNLESS\n")
cat("  H2O's JVM was not a tracked child process (see me_profile.R's header\n")
cat("  caveat and the unfiltered process list it prints) - in that case this\n")
cat("  figure is an UNDERestimate and should be cross-checked hard against\n")
cat("  `qstat -fx <jobid> | grep resources_used` on the first real subjobs,\n")
cat("  not just spot-checked.\n")
cat(sprintf("\n#PBS -l walltime=%s\n", walltime_str))
cat(sprintf("#PBS -l select=1:ncpus=%d:ompthreads=%d:mem=%dgb\n", ncpus, ncpus, mem_gb))
cat(sprintf("Rscript me_analysis.R \"$PBS_ARRAY_INDEX\" %d %d\n\n",
            best$workers, best$n_cores))
cat("Per-task sizing is only half the job: the array's CONCURRENCY throttle\n")
cat("(the %N in #PBS -J 1-360%N) still needs setting from the queue's real\n")
cat("memory/fair-share limits - each concurrent task starts its own H2O JVM\n")
cat("cluster, which this script does not model. See README.md.\n")

# ---- write the directives into me_1.sh --------------------------------------
# workers/n_cores are written as trailing args on the Rscript line, not as a
# manual edit to me_analysis.R's defaults, so the PBS resource request and
# the R-level parallelism config can never drift apart.
if (file.exists(jobscript)) {
  lines <- readLines(jobscript)
  lines <- sub("^#PBS -l walltime=.*$",
               sprintf("#PBS -l walltime=%s", walltime_str), lines)
  lines <- sub("^#PBS -l select=.*$",
               sprintf("#PBS -l select=1:ncpus=%d:ompthreads=%d:mem=%dgb",
                       ncpus, ncpus, mem_gb), lines)
  lines <- sub('^Rscript me_analysis\\.R "\\$PBS_ARRAY_INDEX".*$',
               sprintf('Rscript me_analysis.R "$PBS_ARRAY_INDEX" %d %d',
                       best$workers, best$n_cores),
               lines)
  writeLines(lines, jobscript)
  cat(sprintf("directives and workers/n_cores written into %s\n", jobscript))
} else {
  warning("me_1.sh not found - directives printed only")
}
