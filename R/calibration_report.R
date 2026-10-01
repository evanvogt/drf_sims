##########
# title: bW calibration report
##########
# Prints, for every scenario at every n a study generates at, the calibrated
# bW, the true ATE, and the ATE's planned and realised power, for both outcome
# types - and for binary outcomes the worst-case treated risks. The design is
# in the "Outcome model and bW calibration" sections of sample_size/continuous/README.md
# and sample_size/binary/README.md: each trial is planned for TARGET_POWER under
# homogeneity, and the true ATE equals the planned effect.
#
#   Rscript R/calibration_report.R
#
# Deterministic and quick - no simulation, no RNG - so it runs anywhere,
# including the cluster. Run it there after any change to the scenario tables
# or calibrate_bW(): bW is rounded to 2 dp, and R 4.3.2 there vs 4.5.3 here
# could round a borderline value differently, so the bW tables should match a
# local run exactly.
#
# Continuous power is for an unadjusted two-sample t-test with n / 2 per arm:
#   planned   the no-heterogeneity SD the calibration assumes, so it differs
#             from TARGET_POWER only through bW's rounding
#   realised  adds the treated arm's Var(g), pooled across the two arms, which
#             the plan knows nothing about
#   MNAR-Y0   adds bU^2 sU^2 on top in BOTH arms - the unobserved U's
#             contribution to the control mean (missing-data studies only)
#   MNAR-tau  adds bU^2 sU^2 on top in the treated arm only - U's contribution
#             to the treatment effect (missing-data studies only)
# The missing-data sets' covariates are correlated, so their planned SD also
# carries 2 b1 b2 Cov(X1, X2) (calibrate_bW()).
#
# Binary power is for a two-proportion test with n / 2 per arm, on the true
# marginal risks. Heterogeneity cannot lower it - each arm's variance is fixed
# by its marginal risk - so planned and realised are one column, which differs
# from TARGET_POWER only through bW's rounding. The effect is on the
# risk-difference scale and the MNAR mechanisms' bU * tanh(U) has mean zero, so
# they leave the RD and the power alone; they only widen the risk bounds (the
# treated risk's the same under either).
#
# Binary floor / ceiling: the lowest and highest treated risk over the whole
# covariate support (treated_risk_bounds()). RD_SCALE keeps them inside
# [RD_EPS, 1 - RD_EPS], give or take bW's 3-dp rounding (5e-4);
# sample_size/binary/bin_verify_hte.R checks that, and re-derives RD_SCALE.

suppressPackageStartupMessages(library(here))
suppressMessages(source(here("R", "dgm_scenarios.R")))

MAIN_NS <- c(100, 250, 500, 1000)
MISSING_N <- 500
# validation/continuous/ splits n = 1000 at interim_prop 0.25-0.75
VALIDATION_SCENARIO <- 2
VALIDATION_NS <- seq(250, 750, by = 50)

#' bW, true ATE and t-test powers for one continuous scenario row at one n
continuous_row <- function(p, n) {
  bW <- calibrate_bW(p, n, "t")
  g <- te_moments(p)
  var0 <- p$b1^2 * p$X1_prob * (1 - p$X1_prob) + p$b2^2 * p$s2^2 + p$s_err^2
  if (is_correlated(p)) var0 <- var0 + 2 * p$b1 * p$b2 * cov_X1X2(p)
  ate <- bW + g$mean
  # extra_var is the treated arm's extra variance, pooled over the two arms
  power_at <- function(extra_var) {
    power.t.test(n = n / 2, delta = abs(ate), sd = sqrt(var0 + extra_var / 2))$power
  }
  u_var <- if (is.null(p$bU)) NA_real_ else p$bU^2 * p$sU^2
  c(bW = bW, ate = ate, planned = power_at(0), realised = power_at(g$var),
    mnar_y0 = if (is.na(u_var)) NA_real_ else power_at(g$var + 2 * u_var),
    mnar_tau = if (is.na(u_var)) NA_real_ else power_at(g$var + u_var))
}

#' bW, true RD, two-proportion-test power and worst-case treated risks for one
#' binary scenario row at one n
binary_row <- function(p, n) {
  bW <- calibrate_bW(p, n, "prop")
  p0 <- control_event_rate(p)
  ate <- bW + te_moments(p)$mean
  b <- treated_risk_bounds(p, bW)
  b_y <- if (is.null(p$bU)) c(NA_real_, NA_real_) else treated_risk_bounds(p, bW, p$bU)
  c(bW = bW, p0 = p0, ate = ate,
    power = power.prop.test(n = n / 2, p1 = p0, p2 = p0 + ate)$power,
    floor = b[1], ceiling = b[2],
    floor_mnar = b_y[1], ceiling_mnar = b_y[2])
}

#' One row function over every scenario of a table at every n in MAIN_NS
run_table <- function(tbl, row_fn) {
  width <- length(row_fn(tbl[1, ], MAIN_NS[1]))
  lapply(tbl$scenario, function(s) {
    vapply(MAIN_NS, function(n) row_fn(tbl[tbl$scenario == s, ], n), numeric(width))
  })
}

show <- function(rows, tbl, set, what, label, digits) {
  m <- t(vapply(rows, function(r) r[what, ], numeric(length(MAIN_NS))))
  dimnames(m) <- list(paste("scenario", tbl$scenario), paste0("n=", MAIN_NS))
  cat("\n=== ", set, ": ", label, " ===\n", sep = "")
  print(round(m, digits))
}

#' One table's rows at MISSING_N, scenario 1 without MNAR-tau (not in the grid)
missing_table <- function(tbl, row_fn, cols, labels, bW_digits = 2) {
  m <- t(vapply(tbl$scenario, function(s) {
    row_fn(tbl[tbl$scenario == s, ], MISSING_N)[cols]
  }, numeric(length(cols))))
  dimnames(m) <- list(paste("scenario", tbl$scenario), labels)
  m[, 1] <- round(m[, 1], bW_digits)
  m[, -1] <- round(m[, -1], 3)
  m[1, grepl("MNAR-tau", labels)] <- NA
  m
}

cat(sprintf("Planned power: %.0f%% (TARGET_POWER), n / 2 per arm\n", 100 * TARGET_POWER))

# ---- continuous -------------------------------------------------------------

tbl <- resolve_set("continuous")
rows <- run_table(tbl, continuous_row)
show(rows, tbl, "continuous", "bW", "bW", 2)
show(rows, tbl, "continuous", "ate", "true ATE (bW + E[g])", 3)
show(rows, tbl, "continuous", "planned", "planned power", 3)
show(rows, tbl, "continuous", "realised", "realised power", 3)

cat(sprintf("\n=== continuous scenario %d: bW at the validation study's stage sizes ===\n",
            VALIDATION_SCENARIO))
pv <- tbl[tbl$scenario == VALIDATION_SCENARIO, ]
v <- vapply(VALIDATION_NS, function(n) calibrate_bW(pv, n, "t"), numeric(1))
names(v) <- paste0("n=", VALIDATION_NS)
print(v)

cat(sprintf(paste0("\n=== continuous_missing at n = %d, correlated covariates ",
                   "(scenario 1 has no MNAR-tau) ===\n"), MISSING_N))
print(missing_table(resolve_set("continuous_missing"), continuous_row,
                    c("bW", "ate", "planned", "realised", "mnar_y0", "mnar_tau"),
                    c("bW", "true ATE", "planned", "realised", "MNAR-Y0", "MNAR-tau")))

# ---- binary -----------------------------------------------------------------

tbl <- resolve_set("binary")
rows <- run_table(tbl, binary_row)
cat(sprintf(paste0("\nbinary control risk m0 = %.2f + %.2f * plogis(b0 + b1 X1 + b2 X2), ",
                   "in [%.2f, %.2f]; control event rate E[m0] = %.3f in every scenario\n"),
            tbl$p0_lo[1], tbl$p0_hi[1] - tbl$p0_lo[1], tbl$p0_lo[1], tbl$p0_hi[1],
            rows[[1]]["p0", 1]))
show(rows, tbl, "binary", "bW", "bW", 3)
show(rows, tbl, "binary", "ate", "true ATE (marginal risk difference)", 3)
show(rows, tbl, "binary", "power", "power (planned = realised)", 3)
show(rows, tbl, "binary", "floor", sprintf("treated-risk floor (>= %.2f)", RD_EPS), 3)
show(rows, tbl, "binary", "ceiling", sprintf("treated-risk ceiling (<= %.2f)", 1 - RD_EPS), 3)

cat(sprintf(paste0("\n=== binary_missing at n = %d, correlated covariates; the ",
                   "MNAR treated-risk bounds hold for MNAR-Y0 and MNAR-tau alike ===\n"),
            MISSING_N))
print(missing_table(resolve_set("binary_missing"), binary_row,
                    c("bW", "p0", "ate", "power", "floor_mnar", "ceiling_mnar"),
                    c("bW", "E[m0]", "true RD", "power", "MNAR floor", "MNAR ceiling"),
                    bW_digits = 3))

# ---- sample_size/correlated/ ------------------------------------------------
# scenarios 1-4 on the copula at each rho in CORR_RHOS; rho = 0 should repeat
# the main tables' rows 1-4. At rho = 0.5 binary scenario 3's floor is 0.007,
# below RD_EPS - the documented exception (R/dgm_scenarios.R)

for (rho in CORR_RHOS) {
  set <- corr_set("continuous", rho)
  tbl <- resolve_set(set)
  rows <- run_table(tbl, continuous_row)
  show(rows, tbl, set, "bW", "bW", 2)
  show(rows, tbl, set, "ate", "true ATE (bW + E[g])", 3)
  show(rows, tbl, set, "realised", "realised power", 3)
}
for (rho in CORR_RHOS) {
  set <- corr_set("binary", rho)
  tbl <- resolve_set(set)
  rows <- run_table(tbl, binary_row)
  show(rows, tbl, set, "bW", "bW", 3)
  show(rows, tbl, set, "ate", "true ATE (marginal risk difference)", 3)
  show(rows, tbl, set, "power", "power (planned = realised)", 3)
  show(rows, tbl, set, "floor", sprintf("treated-risk floor (>= %.2f)", RD_EPS), 3)
  show(rows, tbl, set, "ceiling", sprintf("treated-risk ceiling (<= %.2f)", 1 - RD_EPS), 3)
}
