##########
# title: turn the ci_example profiling runs into PBS directives for cts_miss_ci.sh
##########
# Reads the prof_*.RDS files written by cts_miss_ci_profile.R and prints the
# recommended walltime / mem / (workers1, workers2, grf_threads), the
# observed-against-allocated CPU table, a per-scenario data-generation cost
# breakdown, a per-arm bootstrap cost breakdown, and the CI_boot linear-fit
# diagnostics the extrapolation to the real CI_boot = 200 depends on. A
# structural merge of missing/continuous/cts_miss_profile_summary.R (the
# outer/inner flattening and total_elapsed = data_gen_elapsed + elapsed, since
# this study's data generation - mice(m=50) - is likewise
# CI_boot/workers-independent and additive) and
# crossfitting/confidence_intervals/cf_ci_profile_summary.R (the per-configuration
# elapsed ~ CI_boot fit and extrapolation, for the same reason that script has
# one: profiling at the real CI_boot=200 would cost as much as the production
# run it's sizing). See those two files for the parts unchanged here.
#
# Why extrapolate at all: each bootstrap draw (future_map() iteration inside
# cf_oob_half_boot()/rf_oob_half_boot(), R/bootstrap_ci.R) is an independent
# refit on a half-sample, and results are only pooled (combine_mi_ci()) after
# future_map() returns - so elapsed time should scale ~linearly in CI_boot,
# while peak memory (governed by how many draws/imputations are concurrent
# across workers1/workers2, not by CI_boot) should be roughly flat.
# cts_miss_ci_profile.R profiles at CI_boot in {20, 60} rather than the real
# 200 for exactly this reason; this script fits a line through those points
# per (workers1, workers2, grf_threads) configuration and reports it, so a bad
# extrapolation is visible rather than silently trusted.

# libraries
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(paletteer)
library(here)

# paths
path <- here()
prof_path <- file.path(dirname(path), "results", "missing", "ci_example", "profiling")
jobscript <- here("missing", "ci_example", "jobscripts", "cts_miss_ci.sh")

# safety factors
walltime_factor <- 1.75   # a PBS job killed at the limit saves nothing
mem_factor <- 1.5
target_CI_boot <- 200     # cts_miss_ci_analysis.R's / cts_miss_ci.sh's actual production CI_boot

# ---- read the profiling runs -----------------------------------------------
prof_files <- list.files(prof_path, pattern = "^prof_\\d+\\.RDS$", full.names = TRUE)
stopifnot(length(prof_files) > 0)

profs <- map(prof_files, readRDS)

# one row per (index, workers1, workers2, grf_threads, CI_boot) - the outer
# fields are repeated across each outer row's 16 inner cells. Cells where
# cts_miss_ci_profile.R's tryCatch caught a syrup failure (cell$ok == FALSE,
# see that script's header on repeated syrup::syrup() calls per session) are
# dropped here rather than producing NA rows - they carry no usable
# timing/memory data.
skipped_cells <- 0L
runs <- map_dfr(profs, function(p) {
  map_dfr(p$inner, function(cell) {
    if (!isTRUE(cell$ok)) {
      skipped_cells <<- skipped_cells + 1L
      return(NULL)
    }
    tibble(index = p$index, scenario = p$scenario, run = p$run,
           data_gen_elapsed = p$data_gen_elapsed,
           workers1 = cell$workers1, workers2 = cell$workers2,
           grf_threads = cell$grf_threads, CI_boot = cell$CI_boot,
           elapsed = cell$elapsed, cpu_seconds = cell$cpu_seconds,
           peak_rss_gb = cell$peak_rss_gb, n_procs = cell$n_procs,
           expected_procs = cell$expected_procs,
           median_pct_cpu = cell$median_pct_cpu,
           allocated_pct_cpu = cell$allocated_pct_cpu,
           available_cores = p$available_cores)
  })
}) %>%
  mutate(total_elapsed = data_gen_elapsed + elapsed)

if (skipped_cells > 0) {
  warning(sprintf("%d inner cells had no usable measurement (syrup failed on the profiling run) - excluded", skipped_cells))
}

if (n_distinct(runs$CI_boot) < 2) {
  warning("fewer than 2 distinct CI_boot values were profiled - the linear ",
          "extrapolation to CI_boot = ", target_CI_boot, " cannot be fit; ",
          "re-run cts_miss_ci_profile.R across its full grid first.")
}

# a cell profiled on a machine with fewer cores than it asked for is measuring
# that machine's contention, not the configuration
undersized <- runs %>% filter(available_cores < workers1 * workers2 * grf_threads)
if (nrow(undersized) > 0) {
  warning(sprintf("%d profiling cells ran on a machine with fewer cores than the ",
                  nrow(undersized)),
          "configuration requested - their timings are not usable. ",
          "indices: ", paste(unique(undersized$index), collapse = ", "))
  print(as.data.frame(undersized[, c("index", "workers1", "workers2", "grf_threads", "available_cores")]),
        row.names = FALSE)
}

# future falls back to sequential at workers==1 (at either level), so those
# cells have fewer child processes than the fully-nested case
bad_procs <- runs %>% filter(n_procs != expected_procs)
if (nrow(bad_procs) > 0) {
  warning("some cells did not track the expected number of processes (see ",
          "cts_miss_ci_profile.R's expected_procs hypothesis) - the memory ",
          "figures for those may be wrong. indices: ",
          paste(unique(bad_procs$index), collapse = ", "))
}

# ---- cost by scenario (the axis this study's data-gen cost actually varies on) ---
# method is fixed to multiple_imputation throughout cts_miss_ci_config.R's
# grid, so unlike cts_miss_profile_summary.R's "cost by method" table, the
# only thing left to break data-gen cost down by is scenario
cat("\n=== data generation + mice(m=50) imputation cost by scenario (pooled over run) ===\n")
by_scenario <- runs %>%
  distinct(index, scenario, data_gen_elapsed) %>%
  group_by(scenario) %>%
  summarise(n_runs = n(),
            mean_data_gen_s = mean(data_gen_elapsed),
            max_data_gen_s = max(data_gen_elapsed),
            .groups = "drop") %>%
  arrange(desc(mean_data_gen_s))
print(as.data.frame(by_scenario), digits = 3, row.names = FALSE)

cat("\n=== timings by configuration (observed, at each profiled CI_boot; pooled over scenario/run) ===\n")
config_observed <- runs %>%
  group_by(workers1, workers2, grf_threads, CI_boot) %>%
  summarise(n_runs = n(),
            mean_elapsed = mean(elapsed),
            max_elapsed = max(elapsed),
            mean_peak_rss_gb = mean(peak_rss_gb),
            max_peak_rss_gb = max(peak_rss_gb),
            .groups = "drop") %>%
  arrange(workers1, workers2, grf_threads, CI_boot)
print(as.data.frame(config_observed), digits = 3, row.names = FALSE)

# ---- linear fit: elapsed ~ CI_boot, per (workers1, workers2, grf_threads) configuration ----
# pooled over scenario/run, so this is a real fit (n = distinct CI_boot values
# x scenarios x runs), not just a line through 2 means
cat("\n=== elapsed ~ CI_boot fit, per configuration ===\n")
elapsed_fits <- runs %>%
  group_by(workers1, workers2, grf_threads) %>%
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
          "CI_boot value to check. Affected (workers1, workers2, grf_threads): ",
          paste(sprintf("(%d,%d,%d)", poor_fit$workers1, poor_fit$workers2, poor_fit$grf_threads), collapse = ", "))
}

# extrapolated model-fitting+bootstrap elapsed at the real CI_boot, and the
# cpu-seconds that follows from it - this is what the throttled array is
# actually limited by, on top of the (CI_boot-independent) data-gen cost
config <- elapsed_fits %>%
  mutate(mean_elapsed = pmax(intercept + slope * target_CI_boot, 0),
         mean_cpu_seconds = mean_elapsed * workers1 * workers2 * grf_threads) %>%
  left_join(
    runs %>% group_by(workers1, workers2, grf_threads) %>%
      summarise(max_peak_rss_gb = max(peak_rss_gb), .groups = "drop"),
    by = c("workers1", "workers2", "grf_threads")
  ) %>%
  select(workers1, workers2, grf_threads, mean_elapsed, mean_cpu_seconds, max_peak_rss_gb, r_squared) %>%
  arrange(mean_cpu_seconds)

cat(sprintf("\n=== timings extrapolated to CI_boot = %d ===\n", target_CI_boot))
print(as.data.frame(config), digits = 3, row.names = FALSE)

# ---- oversubscription -------------------------------------------------------
cat("\n=== CPU: observed against allocated ===\n")
cpu <- runs %>%
  group_by(workers1, workers2, grf_threads) %>%
  summarise(allocated_pct = first(allocated_pct_cpu),
            observed_pct = median(median_pct_cpu, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(ratio = observed_pct / allocated_pct,
         verdict = case_when(ratio > 1.1 ~ "oversubscribed",
                             ratio < 0.6 ~ "under-used",
                             TRUE ~ "ok")) %>%
  arrange(workers1, workers2, grf_threads)
print(as.data.frame(cpu), digits = 3, row.names = FALSE)

# the array runs throttled (%43, per cts_miss_ci.sh), so it is throughput
# bound: pick the configuration that minimises cpu-seconds charged at the real
# CI_boot, not raw elapsed time
best <- config %>% slice_min(mean_cpu_seconds, n = 1)
cat(sprintf("\nchosen configuration: workers1=%d, workers2=%d, grf_threads=%d (fewest cpu-seconds at CI_boot=%d)\n",
            best$workers1, best$workers2, best$grf_threads, target_CI_boot))

best_runs <- runs %>% filter(workers1 == best$workers1, workers2 == best$workers2, grf_threads == best$grf_threads)

# ---- per-arm cost: fit vs. CI_boot, extrapolated, at the chosen configuration ----
# generic over whatever arm_times$arm keys are actually present
# (causal_forest/dr_random_forest/dr_oracle/dr_semi_oracle) rather than
# assuming fixed names - mirrors cf_ci_profile_summary.R's boot_fits, but
# against the whole per-arm time (nuisance + stage2 + bootstrap all summed
# across 50 imputations, per cts_miss_ci_profile.R's header), since
# verbose_timing's $timings here isn't broken out into a separate boot-only
# component the way this repo's other studies' bootstrap timings are
cat(sprintf("\n=== per-arm time (summed across 50 imputations): fit vs. CI_boot, extrapolated to %d, at the chosen configuration ===\n",
            target_CI_boot))
arm_fits <- map_dfr(profs, function(p) {
  cell <- keep(p$inner, ~ isTRUE(.x$ok) && .x$workers1 == best$workers1 &&
                 .x$workers2 == best$workers2 && .x$grf_threads == best$grf_threads)
  if (length(cell) == 0) return(NULL)
  map_dfr(cell, ~ mutate(.x$arm_times, CI_boot = .x$CI_boot))
}) %>%
  group_by(arm) %>%
  group_modify(~ {
    fit <- lm(time_total ~ CI_boot, data = .x)
    tibble(per_draw_sec = unname(coef(fit)[2]), intercept_sec = unname(coef(fit)[1]),
           r_squared = summary(fit)$r.squared,
           extrapolated_sec = pmax(unname(coef(fit)[1] + coef(fit)[2] * target_CI_boot), 0))
  }) %>%
  ungroup() %>%
  arrange(desc(extrapolated_sec))
print(as.data.frame(arm_fits), digits = 3, row.names = FALSE)
cat(sprintf("sum across arms (the 4 CATE arms are fit sequentially per imputation, so this adds up): %.1f s\n",
            sum(arm_fits$extrapolated_sec)))

# ---- per process behaviour --------------------------------------------------
cat("\n=== per process, at the chosen configuration ===\n")
by_process <- map_dfr(profs, function(p) {
  cell <- keep(p$inner, ~ isTRUE(.x$ok) && .x$workers1 == best$workers1 &&
                 .x$workers2 == best$workers2 && .x$grf_threads == best$grf_threads)
  if (length(cell) == 0) return(NULL)
  map_dfr(cell, function(c) mutate(c$by_process, index = p$index))
}) %>%
  group_by(role) %>%
  summarise(n = n(),
            mean_max_rss_gb = mean(max_rss_gb),
            mean_median_pct_cpu = mean(median_pct_cpu),
            mean_max_pct_cpu = mean(max_pct_cpu),
            .groups = "drop")
print(as.data.frame(by_process), digits = 3, row.names = FALSE)

# ---- timeline plot for one representative cell ------------------------------
# the largest profiled CI_boot at the chosen config (so the bootstrap loop's
# stage is the most visible part of the timeline), on the slowest scenario/run.
# Wrapped in tryCatch, as cts_miss_profile_summary.R does: the recommendation
# and jobscript rewrite below are the actual deliverable of this script, so a
# plotting failure (e.g. an empty usage tibble - the local-Windows syrup gap)
# must not stop the script before it gets there.
tryCatch({
  rep_cell_info <- map_dfr(profs, function(p) {
    cell <- keep(p$inner, ~ isTRUE(.x$ok) && .x$workers1 == best$workers1 &&
                   .x$workers2 == best$workers2 && .x$grf_threads == best$grf_threads)
    if (length(cell) == 0) return(NULL)
    map_dfr(cell, function(c) tibble(index = p$index, CI_boot = c$CI_boot, elapsed = c$elapsed))
  })
  rep_index <- rep_cell_info$index[which.max(rep_cell_info$CI_boot * 1e6 + rep_cell_info$elapsed)]
  rep_prof <- profs[[which(map_dbl(profs, "index") == rep_index)]]
  rep_cell <- keep(rep_prof$inner, ~ isTRUE(.x$ok) && .x$workers1 == best$workers1 &&
                      .x$workers2 == best$workers2 && .x$grf_threads == best$grf_threads &&
                      .x$CI_boot == max(rep_cell_info$CI_boot))[[1]]

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
    labs(title = sprintf("Resource timeline: workers1=%d, workers2=%d, grf_threads=%d, scenario %d, CI_boot=%d",
                         best$workers1, best$workers2, best$grf_threads, rep_prof$scenario, rep_cell$CI_boot),
         subtitle = "One line per process; sampled by syrup. Data generation/imputation (not shown) happens before this window.",
         x = "Seconds since model-fitting started", y = NULL, colour = NULL)
  ggsave("cts_miss_ci_profile_timeline.png", plot = timeline_plot, path = prof_path,
         width = 21, height = 15, units = "cm")
  cat(sprintf("\ntimeline plot written to %s\n", file.path(prof_path, "cts_miss_ci_profile_timeline.png")))
}, error = function(e) {
  warning("timeline plot failed, skipping it - the recommendation below is unaffected: ",
          conditionMessage(e))
})

# ---- the recommendation -----------------------------------------------------
# total_elapsed-based (data-gen + model-fitting/bootstrap), not elapsed alone,
# since the serial mice() step happens inside the same PBS walltime budget.
# The model-fitting/bootstrap part is the extrapolated value at the real
# CI_boot; data-gen is CI_boot-independent so the observed max is used as-is.
max_elapsed <- max(best_runs$data_gen_elapsed) + best$mean_elapsed[1]
peak_gb <- best$max_peak_rss_gb[1]   # observed directly - memory doesn't scale with CI_boot
ncpus <- best$workers1 * best$workers2 * best$grf_threads

# peak_gb can come back non-finite if every cell at the chosen configuration
# had an empty usage tibble (the local-Windows syrup process-tracking gap -
# see memory: syrup_local_flakiness; on a real HPC node with a working ps,
# this shouldn't happen). Refuse to size or write a nonsensical mem= from
# that rather than let sprintf("%d", -Inf) error out below, or silently
# write "mem=-Infgb" into the real jobscript.
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
cat(sprintf("slowest profiled replicate (data-gen + model-fitting/bootstrap, extrapolated to CI_boot=%d): %.1f s\n",
            target_CI_boot, max_elapsed))
cat(sprintf("peak memory (summed RSS across the process tree, model-fitting/bootstrap phase): %.2f gb\n", peak_gb))
cat("  this overcounts shared library pages, so it is an upper bound - the safe\n")
cat("  direction for a mem= request. cross-check once against qstat -fx on the\n")
cat("  first real subjobs; PBS's cgroup figure should sit below it.\n")
cat(sprintf("\n#PBS -l walltime=%s\n", walltime_str))
if (is.na(mem_gb)) {
  cat(sprintf("#PBS -l select=1:ncpus=%d:ompthreads=%d:mem=<no usable measurement>\n", ncpus, ncpus))
} else {
  cat(sprintf("#PBS -l select=1:ncpus=%d:ompthreads=%d:mem=%dgb\n", ncpus, ncpus, mem_gb))
}
cat(sprintf("Rscript cts_miss_ci_analysis.R \"$PBS_ARRAY_INDEX\" %d %d %d\n\n",
            best$workers1, best$workers2, best$grf_threads))

# ---- write the directives into cts_miss_ci.sh --------------------------------
# workers1/workers2/grf_threads are written as trailing args on the Rscript
# line, alongside the #PBS -l lines, not as a manual edit to
# cts_miss_ci_analysis.R's defaults - so the PBS resource request and the
# R-level parallelism config can never drift apart. cts_miss_ci_rerun.sh is
# NOT auto-rewritten here (mirrors cf_ci_profile_summary.R's precedent of
# touching only cf_ci_1.sh, not cf_ci_rerun.sh) - keep it in sync by hand if
# these numbers change.
if (is.na(mem_gb)) {
  cat("jobscript NOT rewritten - no usable memory measurement (see warning above)\n")
} else if (file.exists(jobscript)) {
  lines <- readLines(jobscript)
  lines <- sub("^#PBS -l walltime=.*$",
               sprintf("#PBS -l walltime=%s", walltime_str), lines)
  lines <- sub("^#PBS -l select=.*$",
               sprintf("#PBS -l select=1:ncpus=%d:ompthreads=%d:mem=%dgb",
                       ncpus, ncpus, mem_gb), lines)
  lines <- sub('^Rscript cts_miss_ci_analysis\\.R "\\$PBS_ARRAY_INDEX".*$',
               sprintf('Rscript cts_miss_ci_analysis.R "$PBS_ARRAY_INDEX" %d %d %d',
                       best$workers1, best$workers2, best$grf_threads),
               lines)
  writeLines(lines, jobscript)
  cat(sprintf("directives and workers1/workers2/grf_threads written into %s\n", jobscript))
} else {
  warning("cts_miss_ci.sh not found - directives printed only")
}
