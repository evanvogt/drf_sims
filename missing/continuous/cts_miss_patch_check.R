##########
# title: audit the HTE back-fill - missing covariates, continuous outcome
##########
# As missing/binary/bin_miss_patch_check.R, which carries the full explanation.
#
#   Rscript cts_miss_patch_check.R           # full audit (opens suspect files)
#   Rscript cts_miss_patch_check.R shallow   # timings and .tmp orphans only
#
# Run it on the HPC login node, not here: the results only exist there, so a
# local run overwrites the reports with counts from an empty tree. See the same
# note in bin_miss_patch_check.R and in check_all.R.
#
# This study's patch completed - check_all.R has it at 8,800/8,800 - so it is
# also the control for the checker itself. A correct audit reports 88 `ok`, 11
# `mi_by_design` and an empty failed_patch_ids.txt here. Anything else is a
# false positive in check_patch_failed(), not a problem with these results.

library(here)
source(here("missing/continuous/cts_miss_config.R"))
source(here("R", "patch_hte_tests.R"))

args <- commandArgs(trailingOnly = TRUE)

check_patch_failed(study, deep = !("shallow" %in% args))
