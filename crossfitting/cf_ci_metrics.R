##########
# title: metric definitions for the crossfitting CI pilot
##########
# Definitions only - no side effects, same pattern as cf_metrics.R. Reuses
# cate_metrics()/variant_labels/family_labels from cf_metrics.R and
# interval_metrics() from R/metrics.R rather than redefining any of them.

require(dplyr)
require(here)

source(here("crossfitting", "cf_metrics.R"))
source(here("R", "metrics.R"))

# the 5 arms in scope for the bootstrap CI pilot - subset, not redefine, the
# display vocabulary already defined in cf_metrics.R
ci_variant_levels <- intersect(variant_levels, c("dcf", "scf_scf", "scf_scf_new", "cf_dcf", "cf_scf"))
ci_variant_labels <- variant_labels[ci_variant_levels]

#' Point + interval metrics for every arm of one CI-pilot replicate
#'
#' @param sim_res one per-run object as written by cf_ci_analysis.R
#' @param scenario scenario index
#' @return tibble, one row per arm (training-sample tau only - the bootstrap
#'   bands are not computed for tau_test, see cf_ci_analysis.R)
run_ci_metrics <- function(sim_res, scenario) {
  bind_rows(lapply(names(sim_res$arms), function(nm) {
    a <- sim_res$arms[[nm]]

    bind_cols(
      cate_metrics(a$tau, sim_res$truth_tau, scenario),
      interval_metrics(a$hb_lb, a$hb_ub, sim_res$truth_tau)
    ) %>%
      mutate(
        scenario = scenario,
        run = sim_res$run,
        arm = nm,
        family = a$family,
        variant = a$variant,
        time_nuisance = a$time_nuisance,
        time_stage2 = a$time_stage2,
        time_boot = a$time_boot,
        CI_boot = sim_res$CI_boot,
        CI_sf = sim_res$CI_sf,
        .before = 1
      )
  }))
}
