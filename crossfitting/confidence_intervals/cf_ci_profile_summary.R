##########
# title: turn the CI pilot's profiling runs into PBS directives for cf_ci_1.sh
##########
# Reads the prof_*.RDS files written by cf_ci_profile.R and prints the recommended
# walltime / mem / (workers, grf_threads), the observed-against-allocated CPU table,
# a per-arm bootstrap cost breakdown, and the CI_boot linear-fit diagnostics that the
# extrapolation to the pilot's real CI_boot = 200 depends on. Parallel to
# crossfitting/cf_profile_summary.R, which does the same for cf_1.sh - see that
# file for the parts unchanged here (undersized-machine checks, process-count
# checks, timeline plot).
#
# Why extrapolate at all, rather than reading walltime off the sweep directly: each
# bootstrap draw (future_map() iteration in rf_half_boot/cf_half_boot) is an
# independent refit of V forests on a half-sample, and the n x CI_boot draws matrix
# is only assembled after future_map returns - so elapsed time should scale
# ~linearly in CI_boot, while peak memory (governed by concurrent draws across
# `workers`, not by CI_boot) should be roughly flat. cf_ci_profile.R profiles at
# CI_boot in {20, 60} rather than the real 200 for exactly this reason; this script
# fits a line through those points (pooled over scenario/run, so R^2 is a real
# diagnostic and not a 2-point tautology) and reports it, so a bad extrapolation is
# visible rather than silently trusted.

# libraries
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(paletteer)
library(here)

# paths
path <- here()
prof_path <- file.path(dirname(path), "results", "crossfitting_ci", "profiling")
jobscript <- here("crossfitting", "confidence_intervals", "jobscripts", "cf_ci_1.sh")

# safety factors
walltime_factor <- 1.75   # a PBS job killed at the limit saves nothing
mem_factor <- 1.5
target_CI_boot <- 200     # cf_ci_analysis.R's / cf_ci_1.sh's actual pilot CI_boot

# ---- read the profiling runs -----------------------------------------------
prof_files <- list.files(prof_path, pattern = "^prof_\\d+\\.RDS$", full.names = TRUE)
stopifnot(length(prof_files) > 0)

profs <- map(prof_files, readRDS)

runs <- map_dfr(profs, function(p) {
  tibble(index = p$index, scenario = p$scenario, run = p$run,
         workers = p$workers, grf_threads = p$grf_threads, CI_boot = p$CI_boot,
         elapsed = p$elapsed, cpu_seconds = p$cpu_seconds,
         peak_rss_gb = p$peak_rss_gb, n_procs = p$n_procs,
         median_pct_cpu = p$median_pct_cpu, allocated_pct_cpu = p$allocated_pct_cpu,
         available_cores = p$available_cores, n_arms = p$n_arms)
})

if (n_distinct(runs$CI_boot) < 2) {
  warning("fewer than 2 distinct CI_boot values were profiled - the linear ",
          "extrapolation to CI_boot = ", target_CI_boot, " cannot be fit; ",
          "re-run cf_ci_profile.R across its full grid first.")
}

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

bad_procs <- runs %>% mutate(expected = if_else(workers > 1, workers + 1L, 1L)) %>%
  filter(n_procs != expected)
if (nrow(bad_procs) > 0) {
  warning("some cells did not track the expected number of processes - ",
          "the memory figures for those may be wrong. indices: ",
          paste(bad_procs$index, collapse = ", "))
}

cat("\n=== timings by configuration (observed, at each profiled CI_boot) ===\n")
config_observed <- runs %>%
  group_by(workers, grf_threads, CI_boot) %>%
  summarise(n_runs = n(),
            mean_elapsed = mean(elapsed),
            max_elapsed = max(elapsed),
            mean_peak_rss_gb = mean(peak_rss_gb),
            max_peak_rss_gb = max(peak_rss_gb),
            .groups = "drop") %>%
  arrange(workers, grf_threads, CI_boot)
print(as.data.frame(config_observed), digits = 3, row.names = FALSE)

# ---- linear fit: elapsed ~ CI_boot, per configuration -----------------------
# pooled over scenario/run, so this is a real fit (n = distinct CI_boot values x
# scenarios x runs), not just a line through 2 means
cat("\n=== elapsed ~ CI_boot fit, per configuration ===\n")
elapsed_fits <- runs %>%
  group_by(workers, grf_threads) %>%
  group_modify(~ {
    fit <- lm(elapsed ~ CI_boot, data = .x)
    tibble(intercept = unname(coef(fit)[1]), slope = unname(coef(fit)[2]),
           r_squared = summary(fit)$r.squared, n_obs = nrow(.x))
  }) %>%
  ungroup()
print(as.data.frame(elapsed_fits), digits = 3, row.names = FALSE)

poor_fit <- elapsed_fits %>% filter(r_squared < 0.9)
if (nrow(poor_fit) > 0) {
  warning("elapsed does not scale linearly in CI_boot (R^2 < 0.9) for ",
          nrow(poor_fit), " configuration(s) - the CI_boot = ", target_CI_boot,
          " extrapolation for these is unreliable. Consider profiling a third ",
          "CI_boot value to check. Affected (workers, grf_threads): ",
          paste(sprintf("(%d,%d)", poor_fit$workers, poor_fit$grf_threads), collapse = ", "))
}

# extrapolated elapsed at the pilot's real CI_boot, and the cpu-seconds that
# follows from it - this is what the throttled array is actually limited by
config <- elapsed_fits %>%
  mutate(mean_elapsed = pmax(intercept + slope * target_CI_boot, 0),
         mean_cpu_seconds = mean_elapsed * workers * grf_threads) %>%
  left_join(
    runs %>% group_by(workers, grf_threads) %>%
      summarise(max_peak_rss_gb = max(peak_rss_gb), .groups = "drop"),
    by = c("workers", "grf_threads")
  ) %>%
  select(workers, grf_threads, mean_elapsed, mean_cpu_seconds, max_peak_rss_gb, r_squared) %>%
  arrange(mean_cpu_seconds)

cat(sprintf("\n=== timings extrapolated to CI_boot = %d ===\n", target_CI_boot))
print(as.data.frame(config), digits = 3, row.names = FALSE)

# ---- oversubscription -------------------------------------------------------
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

# the array runs throttled, so it is throughput bound: pick the configuration
# that minimises cpu-seconds charged at the real CI_boot, not raw elapsed time
best <- config %>% slice_min(mean_cpu_seconds, n = 1)
cat(sprintf("\nchosen configuration: workers=%d, grf_threads=%d (fewest cpu-seconds at CI_boot=%d)\n",
            best$workers, best$grf_threads, target_CI_boot))

best_runs <- runs %>% filter(workers == best$workers, grf_threads == best$grf_threads)

# ---- where the time actually goes: nuisance/stage2 vs. bootstrap, per arm ---
cat("\n=== nuisance + stage 2 time by arm, at the chosen configuration (CI_boot-independent) ===\n")
arm_times <- map_dfr(profs, function(p) {
  if (p$workers != best$workers || p$grf_threads != best$grf_threads) return(NULL)
  p$arm_times
}) %>%
  group_by(arm, family, variant) %>%
  summarise(mean_nuisance = mean(time_nuisance),
            mean_stage2 = mean(time_stage2),
            mean_total = mean(time_nuisance + time_stage2),
            .groups = "drop") %>%
  arrange(desc(mean_total))
print(as.data.frame(arm_times), digits = 3, row.names = FALSE)
cat("note: nuisance time is shared between arms that use the same nuisance object,\n")
cat("so summing mean_total over arms double counts - see the replicate elapsed above.\n")

# ---- the actual bottleneck: bootstrap cost per arm, fit and extrapolated ----
cat(sprintf("\n=== bootstrap time by arm: fit vs. CI_boot, extrapolated to %d, at the chosen configuration ===\n",
            target_CI_boot))
boot_fits <- map_dfr(profs, function(p) {
  if (p$workers != best$workers || p$grf_threads != best$grf_threads) return(NULL)
  mutate(p$boot_times, CI_boot = p$CI_boot)
}) %>%
  group_by(arm) %>%
  group_modify(~ {
    fit <- lm(time_boot ~ CI_boot, data = .x)
    tibble(per_draw_sec = unname(coef(fit)[2]), intercept_sec = unname(coef(fit)[1]),
           r_squared = summary(fit)$r.squared,
           extrapolated_sec = pmax(unname(coef(fit)[1] + coef(fit)[2] * target_CI_boot), 0))
  }) %>%
  ungroup() %>%
  arrange(desc(extrapolated_sec))
print(as.data.frame(boot_fits), digits = 3, row.names = FALSE)
cat(sprintf("sum across arms (the bootstrap loop is sequential, so this adds up): %.1f s\n",
            sum(boot_fits$extrapolated_sec)))

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
# the largest profiled CI_boot at the chosen config, so the bootstrap loop's
# stage is the most visible part of the timeline
rep_prof <- profs[[which.max(vapply(profs, function(p) {
  if (p$workers != best$workers || p$grf_threads != best$grf_threads) return(-1)
  p$CI_boot * 1e6 + p$elapsed   # CI_boot first, elapsed as tiebreak
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
  labs(title = sprintf("Resource timeline: workers=%d, grf_threads=%d, scenario %d, CI_boot=%d",
                       rep_prof$workers, rep_prof$grf_threads, rep_prof$scenario, rep_prof$CI_boot),
       subtitle = "One line per process; sampled by syrup",
       x = "Seconds since start", y = NULL, colour = NULL)
ggsave("cf_ci_profile_timeline.png", plot = timeline_plot, path = prof_path,
       width = 21, height = 15, units = "cm")
cat(sprintf("\ntimeline plot written to %s\n", file.path(prof_path, "cf_ci_profile_timeline.png")))

# ---- the recommendation -----------------------------------------------------
max_elapsed <- max(best$mean_elapsed, na.rm = TRUE)  # extrapolated, already at target_CI_boot
peak_gb <- best$max_peak_rss_gb                      # observed directly - memory doesn't scale with CI_boot

walltime_sec <- ceiling(max_elapsed * walltime_factor / 1800) * 1800  # next 30 min
walltime_str <- sprintf("%02d:%02d:00", walltime_sec %/% 3600, (walltime_sec %% 3600) %/% 60)
mem_gb <- ceiling(peak_gb * mem_factor)
ncpus <- best$workers * best$grf_threads

cat("\n=== recommendation ===\n")
cat(sprintf("extrapolated replicate time at CI_boot=%d: %.1f s\n", target_CI_boot, max_elapsed))
cat(sprintf("peak memory (summed RSS across the process tree, observed - CI_boot-independent): %.2f gb\n", peak_gb))
cat("  this overcounts shared library pages, so it is an upper bound - the safe\n")
cat("  direction for a mem= request. cross-check once against qstat -fx on the\n")
cat("  first real subjobs; PBS's cgroup figure should sit below it.\n")
cat(sprintf("\n#PBS -l walltime=%s\n", walltime_str))
cat(sprintf("#PBS -l select=1:ncpus=%d:ompthreads=%d:mem=%dgb\n", ncpus, ncpus, mem_gb))
cat(sprintf("Rscript cf_ci_analysis.R \"$PBS_ARRAY_INDEX\" %d %d %d\n\n",
            target_CI_boot, best$workers, best$grf_threads))

# ---- write the directives into cf_ci_1.sh -----------------------------------
# workers/grf_threads are written as trailing args on the Rscript line, alongside
# CI_boot, not as a manual edit to cf_ci_analysis.R's defaults - so the PBS
# resource request and the R-level parallelism config can never drift apart the
# way cf_ci_1.sh's workers=2 (Rscript line) vs. ncpus=1 (#PBS -l select) had.
if (file.exists(jobscript)) {
  lines <- readLines(jobscript)
  lines <- sub("^#PBS -l walltime=.*$",
               sprintf("#PBS -l walltime=%s", walltime_str), lines)
  lines <- sub("^#PBS -l select=.*$",
               sprintf("#PBS -l select=1:ncpus=%d:ompthreads=%d:mem=%dgb",
                       ncpus, ncpus, mem_gb), lines)
  lines <- sub('^Rscript cf_ci_analysis\\.R "\\$PBS_ARRAY_INDEX".*$',
               sprintf('Rscript cf_ci_analysis.R "$PBS_ARRAY_INDEX" %d %d %d',
                       target_CI_boot, best$workers, best$grf_threads),
               lines)
  writeLines(lines, jobscript)
  cat(sprintf("directives and CI_boot/workers/grf_threads written into %s\n", jobscript))
} else {
  warning("cf_ci_1.sh not found - directives printed only")
}
