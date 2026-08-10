##########
# title: check for failed simulations - model evaluation study
##########
# The grid and the results path come from the study config, so this script
# and the analysis script cannot disagree about what index i means.
# Writes array indices of the missing runs to jobscripts/failed_ids.txt.

library(here)
source(here("model_evaluation/me_config.R"))

check_failed(study)
