##########
# title: unified rerun-campaign status across every simulation study
##########
# Rscript check_all.R
#
# Read-only. Never writes jobscripts/failed_ids.txt and never calls qsub -
# for resubmitting a specific study's failed runs, keep using its own
# <prefix>_check.R -> failed_ids.txt -> jobscripts/<prefix>_rerun.sh loop.
# This script only answers "which study should I look at next" by counting,
# for every study in R/study_registry.R, how many res_sim_*.RDS files exist
# against how many the grid expects.
#
# Results only exist on the HPC (too large to sync here), so this is written
# and syntax-tested locally but run for real on the HPC login node - same
# split as the existing per-study check scripts. See study_registry.R and
# README.md's Status section for what "rerun" means for each study.
#
# Writes check_all_studies.csv and check_all_studies.md next to this script.
# Those two files are committed (not gitignored, unlike failed_ids.txt) so
# `git push` from the HPC + `git pull` here is how progress gets checked
# from this machine without re-running anything.

library(here)
library(dplyr)
source(here("R", "pipeline.R"))
source(here("R", "study_registry.R"))

#' Load a study's `study_config()` object into its own environment, so
#' sourcing 14 config files in one script run can't leak variables between
#' them (some configs define extra top-level constants alongside `study`).
load_study <- function(config_path, config_var) {
  env <- new.env()
  source(here(config_path), local = env)
  if (!exists(config_var, envir = env, inherits = FALSE)) {
    stop("'", config_var, "' not found after sourcing ", config_path)
  }
  get(config_var, envir = env)
}

#' Count res_sim_*.RDS files actually present for a study, without touching
#' failed_ids.txt (that's what check_failed(..., write = TRUE) is for).
count_found <- function(study) {
  cmb <- combos(study)
  sum(vapply(seq_len(nrow(cmb)), function(i) {
    combo <- cmb[i, , drop = FALSE]
    length(list.files(combo_dir(study, combo), pattern = RES_PATTERN))
  }, integer(1)))
}

status_of <- function(found, expected) {
  if (found <= 0) "not_started"
  else if (found >= expected) "complete"
  else "in_progress"
}

check_one <- function(row) {
  if (isTRUE(row$blocked)) {
    return(data.frame(
      study_name = row$study_name, category = row$category,
      expected_jobs = NA_integer_, found_jobs = NA_integer_,
      missing_jobs = NA_integer_, pct_complete = NA_real_,
      status = "blocked", reason = row$reason, stringsAsFactors = FALSE
    ))
  }

  study <- load_study(row$config_path, row$config_var)
  expected <- nrow(study$grid)
  found <- count_found(study)
  missing <- max(expected - found, 0)

  data.frame(
    study_name = row$study_name, category = row$category,
    expected_jobs = expected, found_jobs = found, missing_jobs = missing,
    pct_complete = round(100 * found / expected, 1),
    status = status_of(found, expected), reason = row$reason,
    stringsAsFactors = FALSE
  )
}

results <- bind_rows(lapply(seq_len(nrow(study_registry)), function(i) {
  row <- study_registry[i, ]
  tryCatch(
    check_one(row),
    error = function(e) {
      data.frame(
        study_name = row$study_name, category = row$category,
        expected_jobs = NA_integer_, found_jobs = NA_integer_,
        missing_jobs = NA_integer_, pct_complete = NA_real_,
        status = paste("ERROR:", conditionMessage(e)), reason = row$reason,
        stringsAsFactors = FALSE
      )
    }
  )
})) %>%
  arrange(category, study_name)

# ---- console ----

print(results, row.names = FALSE)

cat("\n=== TOTALS ===\n")
cat("Expected jobs:", sum(results$expected_jobs, na.rm = TRUE), "\n")
cat("Found jobs:   ", sum(results$found_jobs, na.rm = TRUE), "\n")
cat("Missing jobs: ", sum(results$missing_jobs, na.rm = TRUE), "\n")
cat("\nBy status:\n")
print(table(results$status))
cat("\nBy category:\n")
print(table(results$category))

# ---- csv / markdown ----

csv_path <- here("check_all_studies.csv")
md_path  <- here("check_all_studies.md")

write.csv(results, csv_path, row.names = FALSE)

md_table <- function(df) {
  header <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep    <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  rows <- apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(header, sep, rows), collapse = "\n")
}

# Markdown-only ordering: group by status (not_started -> in_progress ->
# complete) so studies still needing work are visually grouped when scanning
# check_all_studies.md, rather than scattered alphabetically by category.
# blocked and the dynamic "ERROR: <msg>" status (unmatched by status_order)
# both sort ahead of everything else, ERROR: first, as the more actionable
# failures. category/study_name stay as the secondary/tertiary sort within
# each status group. Console print() and the CSV keep the original
# category/study_name order (more diff-friendly for git-tracked history).
status_order <- c("blocked", "not_started", "in_progress", "complete")
status_rank <- function(s) {
  rank <- match(s, status_order)
  ifelse(is.na(rank), 0L, rank)
}
results_by_status <- results %>%
  arrange(status_rank(status), category, study_name)

md_lines <- c(
  "# Rerun campaign status",
  "",
  paste0("Last updated: ", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  "",
  md_table(results_by_status),
  "",
  "## Legend",
  "",
  "- **expected_jobs**: n_sims x number of parameter combinations in the study's grid",
  "- **found_jobs**: res_sim_*.RDS files actually present under the study's results directory",
  "- **status**: not_started (0 found), in_progress (0 < found < expected), complete (all found), blocked (currently fails to run, not scanned)",
  "",
  "Generated by `Rscript check_all.R`. For resubmitting a specific study's",
  "missing runs, use its own `<prefix>_check.R` -> `jobscripts/failed_ids.txt`",
  "-> `qsub jobscripts/<prefix>_rerun.sh` loop; this script only reports."
)
writeLines(md_lines, md_path)

cat("\nWrote", csv_path, "and", md_path, "\n")
