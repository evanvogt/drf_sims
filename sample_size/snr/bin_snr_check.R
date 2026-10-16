##########
# title: check for failed simulations - correlated binary SNR check
##########
# Writes array indices of the missing runs to jobscripts/failed_ids.txt and
# points jobscripts/bin_snr_rerun.sh at it.

library(here)
source(here("sample_size", "snr", "bin_snr_config.R"))

check_failed(study)
