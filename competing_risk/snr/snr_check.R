##########
# title: check for failed simulations - competing-risk SNR check
##########
# Writes array indices of the missing runs to jobscripts/failed_ids.txt and
# points jobscripts/snr_rerun.sh at it.

library(here)
source(here("competing_risk", "snr", "snr_config.R"))

check_failed(study)
