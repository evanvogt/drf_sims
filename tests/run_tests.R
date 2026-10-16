##########
# title: quick invariant checks for the shared code
##########
# Rscript tests/run_tests.R        (from anywhere; about a minute)
#
# Base R only, so there is nothing to install. These are not unit tests of every
# function: each one guards a way the shared code can go wrong without anything
# failing - a changed DGM, a grid that no longer matches its jobscript, a
# formula typo in an interval, a metric with the wrong sign. They need no
# results, so they run on a laptop or the login node. The cluster-only
# estimator checks are R/regression_check.R and crossfitting/cf_testing.R.
#
# Exit status is non-zero if any check fails.

suppressPackageStartupMessages(library(here))

failures <- character()
check <- function(name, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) {
    message("    error: ", conditionMessage(e))
    FALSE
  })
  cat(if (ok) "  ok    " else "  FAIL  ", name, "\n", sep = "")
  if (!ok) failures <<- c(failures, name)
}

suppressPackageStartupMessages({
  source(here("R", "utils.R"))
  source(here("R", "pipeline.R"))
  source(here("R", "study_registry.R"))
  source(here("R", "dgm_scenarios.R"))
  source(here("R", "bootstrap_ci.R"))
  source(here("R", "cate_models.R"))
  source(here("R", "metrics.R"))
})

# ---- 1. the DGM delivers the effect it plans for ---------------------------------
# Every simulated trial is planned for 80% power on an ATE, and bW is set so the
# true ATE equals that planned effect in every scenario. The plan is recomputed
# here from first principles and compared with the mean true CATE of a large
# simulated sample, so a change to a table, a coefficient or te_expr that forgets
# to re-derive bW shows up (bug O and bug P were exactly this). Tolerance 0.02:
# the sample ATE of 5e4 units has a Monte Carlo SE under 0.01, and the bugs it
# guards against moved the ATE by 0.1 or more.
cat("
DGM: simulated true ATE equals the planned effect
")
set.seed(1)
for (set in c("continuous", "binary", "continuous_corr_0.5", "binary_corr_0.5")) {
  tbl <- resolve_set(set)
  worst <- 0
  for (k in tbl$scenario) {
    p <- tbl[tbl$scenario == k, ]
    planned <- if (is_binary_set(set)) {
      p0 <- control_event_rate(p)
      power.prop.test(n = 250, p2 = p0, power = TARGET_POWER)$p1 - p0
    } else {
      v <- p$b1^2 * p$X1_prob * (1 - p$X1_prob) + p$b2^2 * p$s2^2 + p$s_err^2 +
        if (is_correlated(p)) 2 * p$b1 * p$b2 * cov_X1X2(p) else 0
      -power.t.test(n = 250, delta = NULL, sd = sqrt(v), power = TARGET_POWER)$delta
    }
    g <- generate_scenario_data(k, 5e4, set, calib_n = 500)
    worst <- max(worst, abs(mean(g$truth$tau) - planned))
  }
  check(sprintf("%s: all %d scenarios within 0.02 of the plan (worst %.4f)",
                set, nrow(tbl), worst), worst < 0.02)
}

# no scenario generates a risk outside [0, 1] at the studies' sample sizes
ok_risk <- TRUE
for (set in c("binary", "binary_corr_0", "binary_corr_0.5")) {
  for (n in c(100, 250, 500, 1000)) for (k in 1:10) {
    ok_risk <- ok_risk && !inherits(try(generate_scenario_data(k, n, set), silent = TRUE),
                                    "try-error")
  }
}
check("binary sets: every scenario at n = 100-1000 generates valid risks", ok_risk)

# hte_mult (the SNR check's lever) and the risk clamp (binary_snr only) leave
# every other study's DGM exactly as it was
same <- TRUE
for (set in c("continuous_corr_0.5", "binary_corr_0.5", "binary_missing")) {
  for (k in 2:4) {
    mech <- if (set == "binary_missing") "MAR" else NULL
    same <- same && identical(
      generate_scenario_data(k, 200, set, mech = mech, seed = k),
      generate_scenario_data(k, 200, set, mech = mech, seed = k, hte_mult = 1)
    ) && identical(get_oracle_info(k, 0.1, set), get_oracle_info(k, 0.1, set, hte_mult = 1))
  }
}
check("hte_mult = 1 is the default DGM, data, truth and oracle alike", same)
check("only binary_snr clamps its risks",
      identical(names(SCENARIO_SETS)[vapply(SCENARIO_SETS, risk_clamped, logical(1))],
                "binary_snr"))

# ---- 2. every first-run jobscript range matches its grid -------------------------
# The array index is a row of study$grid. A -J range that is too short silently
# drops the tail of the study; too long runs jobs with no grid row.
cat("\nJobscripts: -J ranges cover the grid exactly\n")
for (i in seq_len(nrow(study_registry))) {
  reg <- study_registry[i, ]
  s <- load_study(reg$config_path)
  dir <- file.path(dirname(s$failed_file))
  scripts <- list.files(dir, pattern = paste0("^", s$prefix, "_(\\d+|extra)\\.sh$"),
                        full.names = TRUE)
  ranges <- unlist(lapply(scripts, function(f) {
    j <- grep("^#PBS -J ", sub("\r$", "", readLines(f, warn = FALSE)), value = TRUE)
    m <- regmatches(j, regexec("-J (\\d+)-(\\d+)", j))[[1]]
    if (length(m)) seq(as.integer(m[2]), as.integer(m[3])) else integer()
  }))
  check(sprintf("%s: %d jobscript(s) cover rows 1-%d once each",
                reg$study_name, length(scripts), nrow(s$grid)),
        length(scripts) > 0 && identical(sort(ranges), seq_len(nrow(s$grid))))
}

# ---- 3. the simultaneous band has its nominal coverage ----------------------------
# On iid N(0, 1) roots the band should contain a fresh vector of roots about 95%
# of the time. Catches a wrong quantile level (alpha vs alpha/2) or a lost sqrt.
cat("\nBootstrap bands\n")
set.seed(2)
n_u <- 30; B <- 4000
draws <- matrix(rnorm(n_u * B), n_u, B)
band <- simultaneous_band(draws, tau = rep(0, n_u), alpha = 0.05)
fresh <- matrix(rnorm(n_u * 4000), n_u, 4000)
cover <- mean(apply(fresh, 2, function(z) all(z >= band$hb_lb & z <= band$hb_ub)))
check(sprintf("simultaneous_band covers a fresh vector %.3f of the time (want 0.95 +- 0.02)", cover),
      abs(cover - 0.95) < 0.02)

# ---- 4. the three MI strategies use the same margin --------------------------------
# With one imputation the hybrid band must equal simultaneous_band() on its draws.
# It did not until 2026-10-11 (sqrt(lambda * S) instead of sqrt(lambda) * S).
tau <- rnorm(n_u)
one <- combine_mi_ci(list(list(causal_forest = list(tau = tau, draws = draws))),
                     "causal_forest", alpha = 0.05)
ref <- simultaneous_band(draws, tau, alpha = 0.05)
check("combine_mi_ci hybrid band equals simultaneous_band for one imputation",
      isTRUE(all.equal(one$lb_hybrid, ref$hb_lb)) && isTRUE(all.equal(one$ub_hybrid, ref$hb_ub)))

# ---- 5. metric conventions ----------------------------------------------------------
cat("\nMetrics\n")
true <- c(-1, 0, 1, 2)
m <- cate_metrics(true + 0.1, true, scenario = 2)
check("bias and ate_bias are estimate - truth (positive when estimates are too high)",
      isTRUE(all.equal(m$bias, 0.1)) && isTRUE(all.equal(m$ate_bias, 0.1)) &&
        isTRUE(all.equal(m$rmse, 0.1)))
check("indep_p_value: a failed test is NA, not p = 1",
      is.na(indep_p_value(list(p_value = 1, method = "independence_test_failed"))) &&
        identical(indep_p_value(list(p_value = 0.2, method = "independence_test")), 0.2))

# ---- 6. result saving and the stage-2 name clash ---------------------------------------
cat("\nPlumbing\n")
d <- file.path(tempdir(), "save_run_check")
save_run(list(x = 1), d, 7)
save_run(list(x = 2), d, 7)
check("save_run leaves exactly one finished file, overwriting on a repeat",
      identical(list.files(d, all.files = TRUE, no.. = TRUE), "res_sim_7.RDS") &&
        readRDS(file.path(d, "res_sim_7.RDS"))$x == 2)
check("save_run's temporary name is not mistaken for a finished run",
      !grepl(RES_PATTERN, ".tmp_res_sim_7.RDS"))

# sl_plan moves the SuperLearner arm onto its own workers. Its answer must not
# move with it, and the caller's plan must be back afterwards.
sl_fit <- function(sl_plan) {
  setup_rng_stream(3)
  dat <- generate_scenario_data(2, 200, "continuous", seed = 3)$dataset
  suppressWarnings(capture.output(fit <- cate_methods(
    dat, n_folds = 3, sl_lib = c("SL.mean", "SL.glm"), models = "dr_superlearner",
    num.threads = 1, sl_plan = sl_plan)))
  fit
}
plan(sequential)
fit_seq <- sl_fit(NULL)
fit_ms <- sl_fit(tweak(multisession, workers = 2))
check("sl_plan leaves the SuperLearner arm's estimates unchanged and restores the plan",
      length(fit_seq$dr_superlearner$tau) == 200 &&
        identical(fit_seq$dr_superlearner$tau, fit_ms$dr_superlearner$tau) &&
        identical(fit_seq$nuisances_sl, fit_ms$nuisances_sl) &&
        inherits(plan(), "sequential"))

suppressPackageStartupMessages(source(here("crossfitting", "cf_models.R")))
check("sourcing crossfitting/cf_models.R leaves the shared stage2_whole_rf alone",
      all(c("ipw", "Z_query") %in% names(formals(stage2_whole_rf))))

# ---- result ------------------------------------------------------------------------------
cat("\n")
if (length(failures)) {
  cat(length(failures), "check(s) FAILED:\n", paste0("  - ", failures, collapse = "\n"), "\n")
  quit(status = 1)
}
cat("All checks passed.\n")
