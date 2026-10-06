##########
# title: application sizing - which rows still need running
##########
# Writes jobscripts/failed_ids.txt and points as_rerun.sh at it.

library(here)
source(here("application_sizing", "as_config.R"))

check_failed(study)
