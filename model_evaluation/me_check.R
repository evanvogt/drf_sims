##########
# title: check for failed simulations - model evaluation study
##########
# The grid and the results path come from the study config, so this script
# and the analysis script cannot disagree about what index i means.
# Writes array indices of the missing runs to jobscripts/failed_ids*.txt.
#
#   Rscript me_check.R                 # the main study (me_analysis.R)
#   Rscript me_check.R strategies      # the nuisance-arm pass (me_strategies.R)
#   Rscript me_check.R split           # the 80:20 arm (me_split.R)
#
# WHY THE DERIVED TREES ARE CHECKED WITHOUT WRITING. study_strat's grid is all
# 360 rows, because the array index has to keep meaning the same row of the
# same grid as me_analysis.R's. But 2 of those runs have no source file - they
# failed repeatedly in the main study and are deliberately excluded - so
# check_failed() will report them missing forever. Writing that to
# failed_ids_strat.txt would queue them for resubmission, where
# me_strategies.R would find no input and exit 0, and they would be queued
# again next check. So the derived trees are checked with write = FALSE
# (documented in R/pipeline.R as touching nothing) and the success condition
# is "exactly the runs the source tree is also missing", not "none".
#
# The same applies to study_split, whose 240-row grid is the n = 500/1000
# slice: its expected count already excludes n = 250, but not the excluded
# runs that fall inside that slice.

library(here)
source(here("model_evaluation/me_config.R"))

which_study <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(which_study)) which_study <- "main"

st <- me_study(which_study)

if (which_study == "main") {
  check_failed(st)
} else {
  missing_idx <- check_failed(st, write = FALSE)

  # what the source tree is missing, so the two can be compared rather than
  # the derived tree being judged against an unreachable zero
  src_missing <- check_failed(study, write = FALSE)
  src_missing_keys <- do.call(
    paste, c(study$grid[src_missing, c(study$path_cols, "run")], sep = "/")
  )
  missing_keys <- do.call(
    paste, c(st$grid[missing_idx, c(st$path_cols, "run")], sep = "/")
  )

  expected <- nrow(st$grid) - length(intersect(missing_keys, src_missing_keys))
  done <- nrow(st$grid) - length(missing_idx)

  cat(sprintf(
    "\n%s: %d/%d done (%d reachable; %d skipped for want of a source run)\n",
    st$name, done, expected,
    expected, length(intersect(missing_keys, src_missing_keys))
  ))

  unreachable <- setdiff(missing_keys, src_missing_keys)
  if (length(unreachable)) {
    cat("\nmissing WITH a source run available - these are this pass's own",
        "failures and need resubmitting:\n")
    cat(paste0("  ", unreachable, collapse = "\n"), "\n")
    cat("\narray indices:",
        paste(missing_idx[missing_keys %in% unreachable], collapse = " "), "\n")
  } else {
    cat("\nnothing missing that has a source run - this pass is complete.\n")
  }
}
