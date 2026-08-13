##########
# title: back-fill the dr_random_forest HTE tests - missing covariates, binary outcome
##########
# A one-off repair of finished results, not part of the simulation pipeline.
# The work is in R/patch_hte_tests.R; this supplies the study config and the
# array index, exactly as bin_miss_analysis.R does for the simulation itself.
#
#   Rscript bin_miss_patch.R           # every combination, writes
#   Rscript bin_miss_patch.R 7         # combination 7 only (the array element)
#   Rscript bin_miss_patch.R dry       # every combination, writes nothing
#   Rscript bin_miss_patch.R 7 dry     # combination 7, writes nothing
#
# The index is a row of combos(study) - a parameter combination, 100 result
# files - NOT a row of study$grid. There are 99 of them; see
# jobscripts/bin_miss_patch.sh.

library(here)
source(here("missing/binary/bin_miss_config.R"))
source(here("R", "patch_hte_tests.R"))

args <- commandArgs(trailingOnly = TRUE)

# "dry" is matched by value rather than position so the two arguments can be
# given in either order, and so a bare `dry` means "the whole tree, safely".
dry_run <- "dry" %in% args
idx_arg <- setdiff(args, "dry")
combo_idx <- if (length(idx_arg)) as.integer(idx_arg[1]) else NULL

patch_hte_tests(study, combo_idx = combo_idx, dry_run = dry_run)
