##########
# title: check for failed simulations - correlated optimal sample.fraction calibration (continuous)
##########
# The grid and the results path come from the study config, so this script and
# the analysis script cannot disagree about what index i means.
# Writes array indices of the missing runs to jobscripts/failed_cts_ids.txt.

library(here)
source(here("sample_size/correlated/confidence_intervals/optimal_sf/cts_corr_ci_sf_config.R"))

check_failed(study)
