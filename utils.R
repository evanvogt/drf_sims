##########
# title: shim - the shared helpers now live in R/utils.R
##########
# Kept so that existing `source(here("utils.R"))` calls across the studies (and
# in case_study/, validation/, model_evaluation/) keep working unchanged.
# Add new shared helpers to R/utils.R, not here.

source(here::here("R", "utils.R"))
