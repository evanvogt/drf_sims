##########
# title: check for failed simulations - missing covariates, binary outcome
##########
# The grid and the results path come from the study config, so this script and
# the analysis script cannot disagree about what index i means.
# Writes array indices of the missing runs to jobscripts/failed_ids.txt.

library(here)
source(here("missing/binary/bin_miss_config.R"))

check_failed(study)
