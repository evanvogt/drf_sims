##########
# title: optimal sample.fraction calibration for half-sample bootstrap CIs
##########
# find_optimal_sf now lives in R/bootstrap_ci.R alongside the bootstrap it
# calibrates. It is outcome-agnostic despite this file's cts_ prefix - the
# binary calibration study has always sourced it too.

source(here::here("R", "bootstrap_ci.R"))
