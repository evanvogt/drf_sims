##########
# title: turn the profiling runs into PBS directives for cts_val_1.sh
##########
# Reads the prof_*.RDS files written by cts_val_profile.R and prints the
# recommended walltime / mem / (workers, grf_threads), plus the
# observed-against-allocated CPU table, a parallel-efficiency table, a
# cost-against-ncpus table, and a time-by-block table. Everything comes from
# syrup's in-R measurements, so there is no scheduler output to parse. Port of
# continuous/cts_profile_summary.R; see that file for the parts that are
# unchanged.
#
# Run it in the same RStudio session as the sweep:
#
#   source(here::here("validation", "continuous", "cts_val_profile_summary.R"))
#
# Three differences from the continuous version:
#
#   - the cost-against-ncpus table. The point of this exercise is cutting
#     cts_val_1.sh's ncpus=5, so the walltime paid for each core given up has to
#     be visible on its own, not just implied by the (workers, grf_threads)
#     table.
#   - by_interim_prop replaces by_n. n is fixed at 1000 here; what varies is
#     where the trial is split, and the two chunk fits are not symmetric in cost.
#   - the memory line reports peak, session baseline and the difference. The
#     sweep runs every cell in one interactive R session rather than a fresh
#     process per cell, so the raw peak carries whatever that session was
#     already holding. Sizing still uses the raw peak - the safe direction - but
#     a large gap means the recommendation has slack in it.

# libraries
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(paletteer)
library(here)

# paths
source(here("R", "pipeline.R"))        # study_config, read_script, write_script
source(here("validation", "continuous", "cts_val_config.R"))

prof_path <- file.path(study$res_path, "profiling")
jobscript <- here("validation", "continuous", "jobscripts", "cts_val_1.sh")
rerun_script <- here("validation", "continuous", "jobscripts", "cts_val_rerun.sh")

# safety factors
walltime_factor <- 1.75   # a PBS job killed at the limit saves nothing
mem_factor <- 1.5

# ---- read the profiling runs -----------------------------------------------
prof_files <- list.files(prof_path, pattern = "^prof_\\d+\\.RDS$", full.names = TRUE)
stopifnot(length(prof_files) > 0)

profs <- map(prof_files, readRDS)

runs <- map_dfr(profs, function(p) {
  # not `p = p$p`: tibble() evaluates its arguments in sequence against a mask of
  # the columns already built, so a column called `p` shadows the `p` this
  # function was handed and every later p$... reads from an atomic vector
  tibble(index = p$index, interim_prop = p$interim_prop, n1 = p$n1, n2 = p$n2,
         n_covariates = p$p, run = p$run,
         workers = p$workers, grf_threads = p$grf_threads,
         ncpus = p$workers * p$grf_threads,
         elapsed = p$elapsed, cpu_seconds = p$cpu_seconds,
         peak_rss_gb = p$peak_rss_gb, baseline_rss_gb = p$baseline_rss_gb,
         peak_above_baseline_gb = p$peak_above_baseline_gb,
         n_procs = p$n_procs,
         median_pct_cpu = p$median_pct_cpu, allocated_pct_cpu = p$allocated_pct_cpu,
         available_cores = p$available_cores)
})

cat(sprintf("%d profiled cells over %d configurations\n",
            nrow(runs), nrow(distinct(runs, workers, grf_threads))))

# a cell profiled on a machine with fewer cores than it asked for is measuring
# that machine's contention, not the configuration
undersized <- runs %>% filter(available_cores < ncpus)
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

# ---- cost by where the trial is split ---------------------------------------
# n is fixed at 1000, but one job fits two chunks of (n*p, n*(1-p)) and forest
# cost is not linear in n, so the split can matter on its own
cat("\n=== timings by interim_prop (pooled over workers/grf_threads/run) ===\n")
by_interim <- runs %>%
  group_by(interim_prop, n1, n2) %>%
  summarise(n_runs = n(),
            mean_elapsed = mean(elapsed),
            max_elapsed = max(elapsed),
            mean_peak_rss_gb = mean(peak_rss_gb),
            max_peak_rss_gb = max(peak_rss_gb),
            .groups = "drop") %>%
  arrange(interim_prop)
print(as.data.frame(by_interim), digits = 3, row.names = FALSE)

cat("\n=== timings by configuration (pooled over interim_prop/run) ===\n")
config <- runs %>%
  group_by(workers, grf_threads) %>%
  summarise(ncpus = first(ncpus),
            n_runs = n(),
            mean_elapsed = mean(elapsed),
            max_elapsed = max(elapsed),
            mean_cpu_seconds = mean(cpu_seconds),
            mean_peak_rss_gb = mean(peak_rss_gb),
            max_peak_rss_gb = max(peak_rss_gb),
            .groups = "drop") %>%
  arrange(mean_cpu_seconds)
print(as.data.frame(config), digits = 3, row.names = FALSE)

# ---- what each core is actually buying --------------------------------------
# the table this exercise exists for: how much wall time is paid for dropping a
# core, whatever shape the parallelism takes
cat("\n=== cost against ncpus (the core-count question) ===\n")
# at a given core count there can be more than one shape (4 workers x 1 thread
# and 2 x 2 both cost 4 cores), so name the faster one alongside the pooled
# figures rather than hiding the difference in a mean
best_shape <- config %>%
  group_by(ncpus) %>%
  slice_min(mean_elapsed, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(ncpus, best_shape = sprintf("%dw x %dt", workers, grf_threads))

by_ncpus <- runs %>%
  group_by(ncpus) %>%
  summarise(n_runs = n(),
            mean_elapsed = mean(elapsed),
            max_elapsed = max(elapsed),
            mean_cpu_seconds = mean(cpu_seconds),
            max_peak_rss_gb = max(peak_rss_gb),
            .groups = "drop") %>%
  left_join(best_shape, by = "ncpus") %>%
  arrange(ncpus) %>%
  mutate(elapsed_vs_1core = mean_elapsed / mean_elapsed[1]) %>%
  select(ncpus, best_shape, everything())
print(as.data.frame(by_ncpus), digits = 3, row.names = FALSE)
cat(sprintf("\nthe TE-VIM future_map() is over %d covariates, and the only other\n",
            runs$n_covariates[1]))
cat("parallel step is a 9-combination tuning grid - workers beyond that are idle\n")
cat("on those steps and do nothing at all for the sequential spine.\n")

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

# 1100 jobs go out at once with no %N throttle, so this is throughput bound:
# pick the configuration that minimises mean cpu-seconds charged
best <- config %>% slice_min(mean_cpu_seconds, n = 1)
cat(sprintf("\nchosen configuration: workers=%d, grf_threads=%d, ncpus=%d (fewest cpu-seconds)\n",
            best$workers, best$grf_threads, best$ncpus))
cat("  override by hand from the ncpus table above if the wall time it costs is\n")
cat("  not worth the cores it saves.\n")

best_runs <- runs %>% filter(workers == best$workers, grf_threads == best$grf_threads)

# ---- where the time actually goes ------------------------------------------
cat("\n=== mean time by block and chunk, at the chosen configuration ===\n")
# nuisance_rf / causal_forest / dr_random_forest are the sequential spine;
# *_te_vims and *_shap_vims are the only blocks workers touch. That split is
# what decides whether the workers are earning their cores.
arm_times <- map_dfr(profs, function(p) {
  if (p$workers != best$workers || p$grf_threads != best$grf_threads) return(NULL)
  p$arm_times
}) %>%
  group_by(arm, chunk) %>%
  summarise(mean_time = mean(time_total),
            max_time = max(time_total),
            .groups = "drop") %>%
  arrange(desc(mean_time))
print(as.data.frame(arm_times), digits = 3, row.names = FALSE)

cat("\n=== the same, summed over both chunks - one array job's real bill ===\n")
arm_totals <- map_dfr(profs, function(p) {
  if (p$workers != best$workers || p$grf_threads != best$grf_threads) return(NULL)
  p$arm_times %>% mutate(index = p$index)
}) %>%
  group_by(index, arm) %>%
  summarise(both_chunks = sum(time_total), .groups = "drop") %>%
  group_by(arm) %>%
  summarise(mean_time = mean(both_chunks), .groups = "drop") %>%
  mutate(parallel_step = arm %in% c("cf_te_vims", "dr_te_vims",
                                    "cf_shap_vims", "dr_shap_vims"),
         pct = 100 * mean_time / sum(mean_time)) %>%
  arrange(desc(mean_time))
print(as.data.frame(arm_totals), digits = 3, row.names = FALSE)
cat(sprintf("\nparallelisable share of the replicate: %.0f%%\n",
            sum(arm_totals$pct[arm_totals$parallel_step])))
cat("Amdahl's law caps what any number of workers can do against the rest.\n")

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
# the slowest cell at the chosen configuration - the case that should dominate
# the sizing below. Both chunk fits show up in it, back to back.
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
  labs(title = sprintf("Resource timeline: workers=%d, grf_threads=%d, interim_prop=%.2f",
                       rep_prof$workers, rep_prof$grf_threads, rep_prof$interim_prop),
       subtitle = "One line per process; both chunk fits, back to back; sampled by syrup",
       x = "Seconds since start", y = NULL, colour = NULL)
ggsave("cts_val_profile_timeline.png", plot = timeline_plot, path = prof_path,
       width = 21, height = 15, units = "cm")
cat(sprintf("\ntimeline plot written to %s\n",
            file.path(prof_path, "cts_val_profile_timeline.png")))

# ---- the recommendation -----------------------------------------------------
max_elapsed <- max(best_runs$elapsed)
peak_gb <- max(best_runs$peak_rss_gb)
baseline_gb <- max(best_runs$baseline_rss_gb, na.rm = TRUE)

walltime_sec <- ceiling(max_elapsed * walltime_factor / 1800) * 1800  # next 30 min
walltime_str <- sprintf("%02d:%02d:00", walltime_sec %/% 3600, (walltime_sec %% 3600) %/% 60)
mem_gb <- ceiling(peak_gb * mem_factor)
ncpus <- best$ncpus

cat("\n=== recommendation ===\n")
cat(sprintf("slowest profiled replicate (both chunks) : %.1f s\n", max_elapsed))
cat(sprintf("peak memory (summed RSS across the process tree): %.2f gb\n", peak_gb))
cat(sprintf("  of which the RStudio session's own baseline was up to %.2f gb - a\n",
            baseline_gb))
cat("  fresh Rscript under PBS starts lower, so the figure above has slack in it.\n")
cat("  It also overcounts shared library pages, which is the safe direction for a\n")
cat("  mem= request. Cross-check once against qstat -fx on the first real subjobs;\n")
cat("  PBS's cgroup figure should sit below it.\n")
cat(sprintf("\n#PBS -l walltime=%s\n", walltime_str))
cat(sprintf("#PBS -l select=1:ncpus=%d:ompthreads=%d:mem=%dgb\n", ncpus, ncpus, mem_gb))
cat(sprintf("Rscript cts_val_analysis.R \"$PBS_ARRAY_INDEX\" %d %d\n\n",
            best$workers, best$grf_threads))

if (ncpus < 5) {
  cat(sprintf("this drops the request from ncpus=5 to ncpus=%d across 1100 jobs.\n",
              ncpus))
  cat("cts_val_1.sh has no %N throttle on -J, unlike cts_1.sh's %380 - worth\n")
  cat("adding one now that each subjob is cheaper and more of them will run at once.\n\n")
}

# ---- write the directives into the jobscripts -------------------------------
# workers/grf_threads are written as trailing args on the Rscript line, not as a
# manual edit to cts_val_analysis.R's defaults, so the PBS resource request and
# the R-level parallelism config can never drift apart.
#
# read_script()/write_script() rather than readLines()/writeLines(): these
# jobscripts are CRLF, and the readLines idiom the other *_profile_summary.R
# scripts use strips the \r from only the lines it edits, turning a three-line
# change into a whole-file diff. See R/pipeline.R's note.
if (file.exists(jobscript)) {
  sc <- read_script(jobscript)
  sc$lines <- sub("^#PBS -l walltime=.*$",
                  sprintf("#PBS -l walltime=%s", walltime_str), sc$lines)
  sc$lines <- sub("^#PBS -l select=.*$",
                  sprintf("#PBS -l select=1:ncpus=%d:ompthreads=%d:mem=%dgb",
                          ncpus, ncpus, mem_gb), sc$lines)
  sc$lines <- sub('^Rscript cts_val_analysis\\.R "\\$PBS_ARRAY_INDEX".*$',
                  sprintf('Rscript cts_val_analysis.R "$PBS_ARRAY_INDEX" %d %d',
                          best$workers, best$grf_threads),
                  sc$lines)
  write_script(sc, jobscript)
  cat(sprintf("directives and workers/grf_threads written into %s\n",
              basename(jobscript)))
} else {
  warning("cts_val_1.sh not found - directives printed only")
}

# the rerun script's -J and resource lines are left alone: check_failed() ->
# update_rerun_script() recomputes those from cts_val_1.sh on the next
# cts_val_check.R run. Its Rscript line is not recomputed by anything, though,
# and it currently passes no workers/grf_threads at all - so a rerun would
# silently fall back to the defaults in cts_val_analysis.R on whatever
# resources update_rerun_script() had bumped it to.
if (file.exists(rerun_script)) {
  sc <- read_script(rerun_script)
  sc$lines <- sub('^Rscript cts_val_analysis\\.R "\\$jobid".*$',
                  sprintf('Rscript cts_val_analysis.R "$jobid" %d %d',
                          best$workers, best$grf_threads),
                  sc$lines)
  write_script(sc, rerun_script)
  cat(sprintf("workers/grf_threads written into %s (-J and resources left for ",
              basename(rerun_script)))
  cat("check_failed())\n")
} else {
  warning("cts_val_rerun.sh not found")
}
