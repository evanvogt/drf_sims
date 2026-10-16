##########
# title: check or collect one study
##########
# Rscript study.R check <config_path>     missing runs -> study$failed_file, and
#                                         the rerun jobscript's -J and resources
# Rscript study.R collect <config_path>   every run -> <res_path>/<prefix>_all.RDS
#
# <config_path> is relative to the repo root, e.g.
# sample_size/continuous/cts_corr_config.R. model_evaluation/ keeps its own
# me_check.R / me_collect.R, and application_sizing/ and crossfitting/ their
# own *_collect.R.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2 || !args[1] %in% c("check", "collect")) {
  stop("usage: Rscript study.R check|collect <config_path>")
}

library(here)
source(here(args[2]))

if (args[1] == "check") {
  check_failed(study)
} else {
  saveRDS(get_results(study, workers = 2),
          file.path(study$res_path, paste0(study$prefix, "_all.RDS")))
  print("Collection complete!")
}
