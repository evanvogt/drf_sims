##########
# title: check for failed simulations - crossfitting CI pilot
##########
# Mirrors cf_check.R. Kept separate from failed_ids.txt (the production
# study's file) so a re-run of one never clobbers the other's todo list.

library(here)
source(here("crossfitting/confidence_intervals/cf_ci_config.R"))

check_failed(study)
