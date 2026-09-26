##########
# title: continuous bW calibration report
##########
# Prints, for every continuous scenario at every n a study generates at, the
# calibrated bW, the true ATE, and the ATE's planned and realised power. The
# design is in continuous/README.md's "Outcome model and bW calibration": each
# trial is planned for CTS_POWER under homogeneity, the true ATE equals the
# planned effect, and the heterogeneity the plan knows nothing about lowers the
# power actually realised.
#
#   Rscript R/calibration_report.R
#
# Deterministic and quick - no simulation, no RNG - so it runs anywhere,
# including the cluster. Run it there after any change to the scenario tables
# or calibrate_bW(): bW is rounded to 2 dp, and R 4.3.2 there vs 4.5.3 here
# could round a borderline value differently, so the bW tables should match a
# local run exactly.
#
# Power is for an unadjusted two-sample t-test with n / 2 per arm:
#   planned   the no-heterogeneity SD the calibration assumes, so it differs
#             from CTS_POWER only through bW's rounding
#   realised  adds the treated arm's Var(g), pooled across the two arms
#   MNAR-Y    adds bU^2 sU^2 on top - the unobserved U's contribution to the
#             treated arm under that mechanism (missing-data studies only)

suppressPackageStartupMessages(library(here))
suppressMessages(source(here("R", "dgm_scenarios.R")))

MAIN_NS <- c(100, 250, 500, 1000)
MISSING_N <- 500
# validation/continuous/ splits n = 1000 at interim_prop 0.25-0.75
VALIDATION_SCENARIO <- 3
VALIDATION_NS <- seq(250, 750, by = 50)

#' bW, true ATE and powers for one scenario row at one n
calibration_row <- function(p, n) {
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

cat(sprintf("Planned power: %.0f%% (CTS_POWER), unadjusted two-sample t-test, n / 2 per arm\n",
            100 * CTS_POWER))

# ---- main continuous table --------------------------------------------------

tbl <- resolve_set("continuous")
rows <- lapply(tbl$scenario, function(s) {
  vapply(MAIN_NS, function(n) calibration_row(tbl[tbl$scenario == s, ], n),
         numeric(5))
})

show <- function(what, label, digits) {
  m <- t(vapply(rows, function(r) r[what, ], numeric(length(MAIN_NS))))
  dimnames(m) <- list(paste("scenario", tbl$scenario), paste0("n=", MAIN_NS))
  cat("\n=== continuous: ", label, " ===\n", sep = "")
  print(round(m, digits))
}

show("bW", "bW", 2)
show("ate", "true ATE (bW + E[g])", 3)
show("planned", "planned power", 3)
show("realised", "realised power", 3)

cat(sprintf("\n=== continuous scenario %d: bW at the validation study's stage sizes ===\n",
            VALIDATION_SCENARIO))
p3 <- tbl[tbl$scenario == VALIDATION_SCENARIO, ]
v <- vapply(VALIDATION_NS, function(n) calibrate_bW(p3, n, "t"), numeric(1))
names(v) <- paste0("n=", VALIDATION_NS)
print(v)

# ---- missing-data table -----------------------------------------------------

tbl_m <- resolve_set("continuous_missing")
m <- t(vapply(tbl_m$scenario, function(s) {
  calibration_row(tbl_m[tbl_m$scenario == s, ], MISSING_N)
}, numeric(5)))
colnames(m) <- c("bW", "true ATE", "planned", "realised", "MNAR-Y")
rownames(m) <- paste("scenario", tbl_m$scenario)
m[, "bW"] <- round(m[, "bW"], 2)
m[, -1] <- round(m[, -1], 3)
cat(sprintf("\n=== continuous_missing at n = %d (scenario 1 has no MNAR-Y) ===\n",
            MISSING_N))
m[1, "MNAR-Y"] <- NA
print(m)
