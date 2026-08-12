##########
# title: check for failed simulations - optimal sample.fraction (continuous)
##########
# The grid and the results path come from the study config, so this script and
# the analysis script cannot disagree about what index i means.
# Writes array indices of the missing runs to jobscripts/failed_cts_ids.txt.

library(here)
source(here("confidence_intervals/optimal_sf/cts_ci_sf_config.R"))

check_failed(study)
