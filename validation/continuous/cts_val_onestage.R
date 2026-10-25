##########
# title: interim-analysis validation - one-stage full-trial comparator, continuous
##########
# Tests every candidate covariate's interaction with treatment in the whole
# trial of 1000, with a Holm correction, on the trials cts_val_analysis.R split
# - the comparator for the two-stage variable-importance approach. Everything
# but the arm's generator and set is shared with binary/ - see
# validation/val_onestage.R for the design.
#
# Writes <study$res_path>/cts_val_onestage.RDS, list(tests, trials, summary).
#
#   Rscript cts_val_onestage.R                    # cluster: with the pairing check
#   Rscript cts_val_onestage.R nocheck <out.RDS>  # local smoke run, no results files

library(here)

source(here("R", "utils.R"))
source(here("validation", "continuous", "cts_val_dgms.R"))
source(here("validation", "continuous", "cts_val_config.R"))
source(here("validation", "val_onestage.R"))

args <- commandArgs(trailingOnly = TRUE)
run_check <- !(length(args) >= 1 && args[1] == "nocheck")
out_file <- if (length(args) >= 2) args[2] else
  file.path(study$res_path, "cts_val_onestage.RDS")

run_onestage(study,
             generate = generate_continuous_scenario_data,
             set = corr_set("continuous", unique(study$grid$rho)),
             check = run_check,
             out_file = out_file)
