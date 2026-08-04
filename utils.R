##########
# title: shim - the shared helpers now live in R/utils.R
##########
# This file MUST stay at the repo root. validation/cate_validation_cts.R:18
# sources it by absolute path on the cluster:
#
#   source("/rds/.../live/scripts/drf_sims/utils.R")
#
# so removing it breaks that script, and validation/ is outside the scope of the
# restructure. case_study/ also calls setup_rng_stream() and
# collate_predictions() without sourcing anything itself, relying on whatever
# drives it to have loaded them.
#
# (model_evaluation/ does NOT depend on this file - it sources its own
# sim_utils.R / model_utils.R / nuisance_utils.R / metric_utils.R.)
#
# Every in-repo script now sources R/utils.R directly, so this shim exists only
# for the out-of-scope folders above. Nothing in the repo should acquire a new
# dependency on it.
#
# Add new shared helpers to R/utils.R, not here.

source(here::here("R", "utils.R"))
