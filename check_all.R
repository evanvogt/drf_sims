##########
# title: expected, found and missing runs for every study
##########
# Rscript check_all.R
#
# Read-only: for every study in R/study_registry.R, counts the res_sim_*.RDS
# files present against the rows of its grid, and prints the table. It never
# writes failed_ids.txt or calls qsub - a study's own <prefix>_check.R does
# that. Run it on the HPC login node, where the results are.

library(here)
library(dplyr)
source(here("R", "pipeline.R"))
source(here("R", "study_registry.R"))

count_found <- function(study) {
  cmb <- combos(study)
  sum(vapply(seq_len(nrow(cmb)), function(i) {
    length(list.files(combo_dir(study, cmb[i, , drop = FALSE]), pattern = RES_PATTERN))
  }, integer(1)))
}

check_one <- function(row) {
  study <- load_study(row$config_path)
  expected <- nrow(study$grid)
  found <- count_found(study)
  data.frame(study_name = row$study_name, expected = expected, found = found,
             missing = max(expected - found, 0),
             pct = round(100 * found / expected, 1), error = NA_character_)
}

results <- bind_rows(lapply(seq_len(nrow(study_registry)), function(i) {
  row <- study_registry[i, ]
  tryCatch(check_one(row), error = function(e) {
    data.frame(study_name = row$study_name, expected = NA_integer_,
               found = NA_integer_, missing = NA_integer_, pct = NA_real_,
               error = conditionMessage(e))
  })
}))
if (all(is.na(results$error))) results$error <- NULL

print(results, row.names = FALSE)
cat("\nExpected:", sum(results$expected, na.rm = TRUE),
    " found:", sum(results$found, na.rm = TRUE),
    " missing:", sum(results$missing, na.rm = TRUE), "\n")
