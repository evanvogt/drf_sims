##########
# title: application sizing - collect the timings
##########
# Writes, into study$res_path:
#   as_fits.csv  one row per timed fit (seconds)
#   as_jobs.csv  one row per grid row: n, p, peak memory
# Peak memory comes two ways: the R processes' own high-water marks (main +
# each future worker, summed - an upper bound, as the peaks need not
# coincide), and, where the PBS log is still in jobscripts/logs*/, the
# scheduler's own peak for the whole job.

library(here)
source(here("application_sizing", "as_config.R"))

files <- vapply(seq_len(nrow(study$grid)), function(i) {
  file.path(combo_dir(study, study$grid[i, ]), paste0("res_sim_", study$grid$run[i], ".RDS"))
}, character(1))
have <- file.exists(files)
cat(sum(have), "of", length(files), "rows have results\n")

res <- lapply(files[have], readRDS)

fits <- do.call(rbind, lapply(res, `[[`, "fits"))

jobs <- do.call(rbind, lapply(which(have), function(i) {
  r <- res[[match(i, which(have))]]
  data.frame(row = i, study$grid[i, c("dataset", "task", "handling", "sf")],
             n = r$n, p = r$p, seconds = sum(r$fits$seconds),
             peak_r_mb = r$peak_main_mb + sum(r$peak_workers_mb, na.rm = TRUE),
             node = r$node, row.names = NULL)
}))

# PBS peak memory, if the logs are here. The job summary footer reads
# "Used : <mem> (peak)"; logs are <jobname>.o<jobid>.<array index>. A rerun
# log's index is a line of failed_ids.txt, not a grid row, so rerun logs are
# skipped - their rows keep the R-side estimate only.
log_files <- list.files(here("application_sizing", "jobscripts"),
                        pattern = "^as_1\\.o", recursive = TRUE, full.names = TRUE)
pbs <- do.call(rbind, lapply(log_files, function(f) {
  txt <- readLines(f, warn = FALSE)
  peak <- regmatches(txt, regexpr("[0-9.]+\\s*\\(peak\\)", txt))
  if (!length(peak)) return(NULL)
  data.frame(row = as.integer(sub(".*\\.", "", basename(f))),
             pbs_peak_gb = as.numeric(sub("\\s*\\(peak\\)", "", peak[1])))
}))
if (!is.null(pbs)) jobs <- merge(jobs, pbs, by = "row", all.x = TRUE)

write.csv(fits, file.path(study$res_path, "as_fits.csv"), row.names = FALSE)
write.csv(jobs, file.path(study$res_path, "as_jobs.csv"), row.names = FALSE)
print("Collection complete!")
