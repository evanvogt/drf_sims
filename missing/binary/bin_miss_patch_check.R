##########
# title: audit the HTE back-fill - missing covariates, binary outcome
##########
# The counterpart to bin_miss_check.R, for the repair rather than the
# simulation. bin_miss_check.R asks "which res_sim_*.RDS are missing"; this asks
# "which of them did bin_miss_patch.sh fail to patch, and why".
#
#   Rscript bin_miss_patch_check.R           # full audit (opens suspect files)
#   Rscript bin_miss_patch_check.R shallow   # timings and .tmp orphans only
#
# WHY IT EXISTS
#
# check_all.R had this study at patchable 8,800 / patched 7,800 - a shortfall of
# exactly 1,000 files. It counts manifest ROWS, and ten of the 99 array elements
# left no manifest in ../results/missing/binary/bin_miss_hte_patch/ (combos 21,
# 22, 28, 29, 30, 32, 40, 44, 45, 48). A missing manifest and a combination that
# was never submitted look identical from there.
#
# The usual next step - read the job's stderr - was not available: the HPC
# returns no .e files and the .o files from that run are gone. So the audit is
# reconstructed from the result files, which do survive; see the "auditing the
# repair" section of R/patch_hte_tests.R for what each trace proves.
#
# The `deep` pass opens every result file in a suspect combination (~100
# readRDS each, a couple of minutes for ten) and is the point of the script:
# only it can tell "the element died before writing its manifest but the files
# ARE patched" from "the files are genuinely unpatched". The first needs a
# re-run purely to regenerate the manifest; the second needs the work redoing.
# `shallow` skips it when you only want the timings.
#
# Writes bin_miss_patch_check.{csv,md} here and jobscripts/failed_patch_ids.txt.
# The reports are committed, unlike results and logs, so `git push` from the HPC
# and `git pull` at the other end is how the answer gets off the cluster.
#
# Which means, exactly as for check_all.R: running this on the laptop OVERWRITES
# all three with counts taken from the two smoke-test result files that live
# there, and reports the study as barely patched. It is a syntax check, not a
# status check - `git checkout --` the three afterwards, or the next commit
# carries the wrong answer to the cluster.

library(here)
source(here("missing/binary/bin_miss_config.R"))
source(here("R", "patch_hte_tests.R"))

args <- commandArgs(trailingOnly = TRUE)

check_patch_failed(study, deep = !("shallow" %in% args))
