##########
# title: shim - the shared helpers now live in R/utils.R
##########
# This file MUST stay at the repo root. case_study/ calls setup_rng_stream()
# and collate_predictions() without sourcing anything itself, relying on
# whatever drives it to have loaded them from here first - case_study/ is
# outside the scope of the restructure.

#
# Every in-repo script now sources R/utils.R directly, so this shim exists only
# for case_study/. Nothing in the repo should acquire a new dependency on it.
#
# Add new shared helpers to R/utils.R, not here.

source(here::here("R", "utils.R"))
