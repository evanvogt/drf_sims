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

# the CI pilot now covers every non-SuperLearner arm, which is exactly what
# cf_metrics.R's variant_levels already lists, in display order. Kept as named
# variables so cf_ci_results.R needs no rename.
ci_variant_levels <- variant_levels
ci_variant_labels <- variant_labels[ci_variant_levels]

# display order and labels for the interval methods. half_boot is available for
# every arm; the other two only for the whole-sample/OOB arms.
ci_method_levels <- c("half_boot", "half_boot_out", "grf_normal")

ci_method_labels <- c(
  half_boot     = "Half-sample bootstrap",
  half_boot_out = "Half-sample bootstrap (out-of-half)",
  grf_normal    = "grf variance (pointwise)"
)

#' Point + interval metrics for every (arm, interval method) of one CI-pilot replicate
#'
#' Several rows per arm is the shape R/metrics.R's compute_metrics already
#' documents for the CI studies, and the reason the extra methods are worth
#' carrying is that they are nearly free: half_boot and half_boot_out come off
#' one set of forest refits, and grf_normal off a variance grf computed while
#' fitting. Only the whole-sample/OOB arms have the latter two.
#'
#' Read simultaneous_coverage with the method in mind: both bootstrap bands are
#' constructed to control it, grf_normal is a pointwise interval and is expected
#' to score near zero on it. That is the contrast, not a defect.
#'
#' @param sim_res one per-run object as written by cf_ci_analysis.R
#' @param scenario scenario index
#' @return tibble, one row per (arm, ci_method) (training-sample tau only - no
#'   interval is computed for tau_test, see cf_ci_analysis.R)
run_ci_metrics <- function(sim_res, scenario) {
  bind_rows(lapply(names(sim_res$arms), function(nm) {
    a <- sim_res$arms[[nm]]
    point <- cate_metrics(a$tau, sim_res$truth_tau, scenario)

    # (method, lb, ub, time_boot) for every interval this arm carries
    intervals <- list(
      list("half_boot", a$hb_lb, a$hb_ub, a$time_boot)
    )
    if (!is.null(a$hb_out_lb)) {
      # same refits as half_boot, so it is charged the same bootstrap time
      intervals <- c(intervals, list(list("half_boot_out", a$hb_out_lb, a$hb_out_ub, a$time_boot)))
    }
    if (!is.null(a$var_oob)) {
      ni <- normal_interval(a$tau, a$var_oob, sim_res$alpha)
      intervals <- c(intervals, list(list("grf_normal", ni$lb, ni$ub, 0)))
    }

    bind_rows(lapply(intervals, function(iv) {
      bind_cols(point, interval_metrics(iv[[2]], iv[[3]], sim_res$truth_tau)) %>%
        mutate(
          scenario = scenario,
          run = sim_res$run,
          arm = nm,
          family = a$family,
          variant = a$variant,
          ci_method = iv[[1]],
          time_nuisance = a$time_nuisance,
          time_stage2 = a$time_stage2,
          time_boot = iv[[4]],
          CI_boot = sim_res$CI_boot,
          CI_sf = sim_res$CI_sf,
          .before = 1
        )
    }))
  }))
}
