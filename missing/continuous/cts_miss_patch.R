##########
# title: back-fill the dr_random_forest HTE tests - missing covariates, continuous outcome
##########
# As missing/binary/bin_miss_patch.R, against this study's config. The work is
# in R/patch_hte_tests.R.
#
#   Rscript cts_miss_patch.R           # every combination, writes
#   Rscript cts_miss_patch.R 7         # combination 7 only (the array element)
#   Rscript cts_miss_patch.R dry       # every combination, writes nothing
#   Rscript cts_miss_patch.R 7 dry     # combination 7, writes nothing
#
# The index is a row of combos(study) - a parameter combination, 100 result
# files - NOT a row of study$grid. There are 99 of them; see
# jobscripts/cts_miss_patch.sh.
#
# Run this AFTER cts_miss_rerun.sh has cleared jobscripts/failed_ids.txt and
# cts_miss_check.R reports 9,900/9,900, so the patch is one clean pass over a
# complete study. It is idempotent, so patching early and again afterwards also
# works - it just splits the manifest across two passes.

library(here)
source(here("missing/continuous/cts_miss_config.R"))
source(here("R", "patch_hte_tests.R"))

args <- commandArgs(trailingOnly = TRUE)

dry_run <- "dry" %in% args
idx_arg <- setdiff(args, "dry")
combo_idx <- if (length(idx_arg)) as.integer(idx_arg[1]) else NULL

patch_hte_tests(study, combo_idx = combo_idx, dry_run = dry_run)
