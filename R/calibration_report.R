##########
# title: bW calibration report
##########
# Prints, for every scenario at every n a study generates at, the calibrated
# bW, the true ATE, and the ATE's planned and realised power, for both outcome
# types. The design is in the "Outcome model and bW calibration" sections of
# continuous/README.md and binary/README.md: each trial is planned for
# TARGET_POWER under homogeneity, and the true ATE equals the planned effect.
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
#   MNAR-Y    adds bU^2 sU^2 on top - the unobserved U's contribution to the
#             treated arm under that mechanism (missing-data studies only)
#
# Binary power is for a two-proportion test with n / 2 per arm, on the true
# marginal risks. Heterogeneity cannot lower it - each arm's variance is fixed
# by its marginal risk - so planned and realised are one column, which differs
# from TARGET_POWER only through bW's rounding. Under MNAR-Y the treated risk
# is also averaged over U, which pulls it towards 0.5 and shrinks the RD.

suppressPackageStartupMessages(library(here))
suppressMessages(source(here("R", "dgm_scenarios.R")))

MAIN_NS <- c(100, 250, 500, 1000)
MISSING_N <- 500
# validation/continuous/ splits n = 1000 at interim_prop 0.25-0.75
VALIDATION_SCENARIO <- 3
VALIDATION_NS <- seq(250, 750, by = 50)

#' bW, true ATE and t-test powers for one continuous scenario row at one n
continuous_row <- function(p, n) {
  bW <- calibrate_bW(p, n, "t")
  g <- te_moments(p)
  var0 <- p$b1^2 * p$X1_prob * (1 - p$X1_prob) + p$b2^2 * p$s2^2 + p$s_err^2
  ate <- bW + g$mean
  power_at <- function(extra_var) {
    power.t.test(n = n / 2, delta = abs(ate), sd = sqrt(var0 + extra_var / 2))$power
  }
  u_var <- if (is.null(p$bU)) NA_real_ else p$bU^2 * p$sU^2
  c(bW = bW, ate = ate, planned = power_at(0), realised = power_at(g$var),
    mnar_y = if (is.na(u_var)) NA_real_ else power_at(g$var + u_var))
}

#' bW, true RD and two-proportion-test powers for one binary scenario row at one n
binary_row <- function(p, n) {
  bW <- calibrate_bW(p, n, "prop")
  base <- baseline_grid(p)
  tg <- te_grid(p)
  p0 <- marginal_risk(base, list(g = 0, w = 1))
  rd <- marginal_risk(base, tg, bW) - p0
  rd_y <- NA_real_
  if (!is.null(p$bU)) {
    # bU * U enters the treated arm's linear predictor alongside b2 * X2, both
    # normal and independent, so averaging over U just widens that term
    wide <- p
    wide$s2 <- sqrt(p$s2^2 + (p$bU * p$sU / p$b2)^2)
    rd_y <- marginal_risk(baseline_grid(wide), tg, bW) - p0
  }
  power_at <- function(d) power.prop.test(n = n / 2, p1 = p0, p2 = p0 + d)$power
  c(bW = bW, p0 = p0, ate = rd, power = power_at(rd),
    rd_mnar_y = rd_y, mnar_y = if (is.na(rd_y)) NA_real_ else power_at(rd_y))
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

#' One table's rows at MISSING_N, scenario 1 without MNAR-Y (not in the grid)
missing_table <- function(tbl, row_fn, cols, labels) {
  m <- t(vapply(tbl$scenario, function(s) {
    row_fn(tbl[tbl$scenario == s, ], MISSING_N)[cols]
  }, numeric(length(cols))))
  dimnames(m) <- list(paste("scenario", tbl$scenario), labels)
  m[, 1] <- round(m[, 1], 2)
  m[, -1] <- round(m[, -1], 3)
  m[1, grepl("MNAR-Y", labels)] <- NA
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
p3 <- tbl[tbl$scenario == VALIDATION_SCENARIO, ]
v <- vapply(VALIDATION_NS, function(n) calibrate_bW(p3, n, "t"), numeric(1))
names(v) <- paste0("n=", VALIDATION_NS)
print(v)

cat(sprintf("\n=== continuous_missing at n = %d (scenario 1 has no MNAR-Y) ===\n",
            MISSING_N))
print(missing_table(resolve_set("continuous_missing"), continuous_row,
                    c("bW", "ate", "planned", "realised", "mnar_y"),
                    c("bW", "true ATE", "planned", "realised", "MNAR-Y")))

# ---- binary -----------------------------------------------------------------

tbl <- resolve_set("binary")
rows <- run_table(tbl, binary_row)
cat(sprintf("\nbinary control-arm risk p0 = E[plogis(b0 + b1 X1 + b2 X2)] = %.3f in every scenario\n",
            rows[[1]]["p0", 1]))
show(rows, tbl, "binary", "bW", "bW", 2)
show(rows, tbl, "binary", "ate", "true ATE (marginal risk difference)", 3)
show(rows, tbl, "binary", "power", "power (planned = realised)", 3)

cat(sprintf("\n=== binary_missing at n = %d (scenario 1 has no MNAR-Y) ===\n",
            MISSING_N))
print(missing_table(resolve_set("binary_missing"), binary_row,
                    c("bW", "ate", "power", "rd_mnar_y", "mnar_y"),
                    c("bW", "true RD", "power", "MNAR-Y RD", "MNAR-Y power")))
