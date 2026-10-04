##########
# title: check for failed simulations - single-event survival
##########
# Writes array indices of the missing runs to jobscripts/failed_ids.txt and
# points se_rerun.sh at them. Run on the cluster, not locally: against a local
# results tree it would rewrite both for every run you have not produced.

library(here)
source(here("competing_risk", "single_event", "se_config.R"))

check_failed(study)
