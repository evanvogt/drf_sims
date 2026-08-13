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
# It also reports any one-off repair a study owes on top of being run, from the
# registry's patch_manifest column. Being run and being correct are different
# states: missing/binary sat at 9,900/9,900 "complete" while still owing the
# dr_random_forest HTE back-fill, and a file count cannot see that.
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
#
# Which means: running this locally OVERWRITES both with local counts, and the
# local machine has almost no results. It is a syntax check, not a status
# check - `git checkout -- check_all_studies.csv check_all_studies.md`
# afterwards, or the next commit reports the whole campaign as barely started.

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

#' How far through its one-off repair a study is
#'
#' Counted from the manifest R/patch_hte_tests.R writes (one small CSV per
#' parameter combination), NOT by opening the result files - reading 9,900 RDS
#' objects to look for one field would take about half an hour per study, which
#' is far too slow for a script meant to be run casually on a login node.
#'
#' Returns patchable (the denominator), patched, and NA/NA for a study that owes
#' no repair. The denominator excludes the multiple_imputation runs: those carry
#' no saved nuisances, so there is nothing to compute a BLP from and the patch
#' rightly refuses them (see the multiple-imputation note in missing/README.md).
#' Counting them would leave a fully repaired study stuck at 88.9% for ever.
count_patched <- function(study, manifest_dir) {
  if (is.na(manifest_dir)) return(list(patchable = NA_integer_, patched = NA_integer_))

  mi_rows <- if ("method" %in% names(study$grid)) {
    length(grid_indices(study, method = "multiple_imputation"))
  } else 0L
  patchable <- nrow(study$grid) - mi_rows

  files <- list.files(file.path(study$res_path, manifest_dir),
                      pattern = "\\.csv$", full.names = TRUE)
  if (length(files) == 0) return(list(patchable = patchable, patched = 0L))

  rows <- bind_rows(lapply(files, read.csv, stringsAsFactors = FALSE))
  # "already_patched" counts: a second pass over a file that was done in the
  # first is still a patched file, and the manifest is per-pass.
  done <- sum(rows$status %in% c("patched", "already_patched"))

  list(patchable = patchable, patched = done)
}

check_one <- function(row) {
  if (isTRUE(row$blocked)) {
    return(data.frame(
      study_name = row$study_name, category = row$category,
      expected_jobs = NA_integer_, found_jobs = NA_integer_,
      missing_jobs = NA_integer_, pct_complete = NA_real_,
      status = "blocked",
      patchable_jobs = NA_integer_, patched_jobs = NA_integer_,
      patch_status = "not_applicable",
      reason = row$reason, stringsAsFactors = FALSE
    ))
  }

  study <- load_study(row$config_path, row$config_var)
  expected <- nrow(study$grid)
  found <- count_found(study)
  missing <- max(expected - found, 0)

  patch <- count_patched(study, row$patch_manifest)

  data.frame(
    study_name = row$study_name, category = row$category,
    expected_jobs = expected, found_jobs = found, missing_jobs = missing,
    pct_complete = round(100 * found / expected, 1),
    status = status_of(found, expected),
    patchable_jobs = patch$patchable, patched_jobs = patch$patched,
    patch_status = if (is.na(patch$patchable)) "not_applicable" else {
      status_of(patch$patched, patch$patchable)
    },
    reason = row$reason,
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
        status = paste("ERROR:", conditionMessage(e)),
        patchable_jobs = NA_integer_, patched_jobs = NA_integer_,
        patch_status = "not_applicable",
        reason = row$reason,
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

# Reported separately from the run counts because a study can be 100% run and
# still owe its repair - which is the whole reason these columns exist.
owed <- results[results$patch_status != "not_applicable", ]
if (nrow(owed)) {
  cat("\nOne-off repairs (R/patch_hte_tests.R):\n")
  print(owed[, c("study_name", "patchable_jobs", "patched_jobs", "patch_status")],
        row.names = FALSE)
}

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
  "- **patchable_jobs / patched_jobs / patch_status**: progress of a one-off",
  "  repair the study owes on top of being run, counted from the manifest the",
  "  repair writes. `not_applicable` means the study owes none. Right now the",
  "  only repair is the dr_random_forest HTE back-fill (`R/patch_hte_tests.R`),",
  "  owed by the two missing-covariate studies. **A study can be `complete` and",
  "  its patch `not_started`** - that is exactly the state these columns exist",
  "  to make visible.",
  "- **patchable_jobs excludes the `multiple_imputation` runs** (1,100 per",
  "  study). Those keep no nuisances, so there is nothing to recompute a BLP",
  "  from and the patch refuses them by design; counting them would peg a fully",
  "  repaired study at 88.9%. See the multiple-imputation note in",
  "  `missing/README.md`.",
  "",
  "Generated by `Rscript check_all.R`. For resubmitting a specific study's",
  "missing runs, use its own `<prefix>_check.R` -> `jobscripts/failed_ids.txt`",
  "-> `qsub jobscripts/<prefix>_rerun.sh` loop; this script only reports."
)
writeLines(md_lines, md_path)

cat("\nWrote", csv_path, "and", md_path, "\n")
