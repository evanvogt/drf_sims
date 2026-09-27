##########
# title: binary study - check the risk-difference scenarios are what they claim
##########
# Every binary set generates P(Y = 1 | x, W) = m0(x) + W * tau(x): the treatment
# effect on the RISK-DIFFERENCE scale, with the control risk m0 bounded in
# [p0_lo, p0_hi] (R/dgm_scenarios.R; the design is in README.md). This script
# checks the properties that design exists for, and re-derives RD_SCALE:
#   1. the oracle formula and the generator's truth agree (they are separate
#      columns of the scenario table and could drift)
#   2. tau depends on the described modifiers only: moving X1 and X2 leaves it
#      unchanged, so scenario 1 is an exact null and the HTE is what te_expr says
#   3. SD(tau) is the same at every n - only the ATE moves - and the ATE and
#      power are the planned ones
#   4. every treated risk stays inside [RD_EPS, 1 - RD_EPS] at every n, and
#      under the missing study's MNAR-Y
#   5. each RD_SCALE[k] is the largest scale check 4 allows, floored to 3 dp
#
# Until 2026-09-26 the effect was on the logit scale, and this script measured
# how much of var(tau) the link handed to X1 and X2 - up to 80% at n = 100.
# README.md keeps that table as the reason for the change.
#
# Writes nothing, and stops at the first failed check. Run from sample_size/binary/:
#   Rscript bin_verify_hte.R

source(here::here("R", "dgm_scenarios.R"))

tbl <- resolve_set("binary")
miss <- resolve_set("binary_missing")
ns <- c(100, 250, 500, 1000)   # binary/; the CI studies run 500 and 1000
MISSING_N <- 500               # the only n in missing/binary/bin_miss_config.R
ROUNDING <- 5e-4               # bW is rounded to 3 dp

check <- function(ok, what) {
  if (!all(ok)) stop("FAILED: ", what, call. = FALSE)
  cat("ok:", what, "\n")
}

# the planned risk difference at n, before bW's rounding - as calibrate_bW()
planned_rd <- function(p, n) {
  p0 <- control_event_rate(p)
  power.prop.test(n = n / 2, p2 = p0, power = TARGET_POWER)$p1 - p0
}

# ---- 1. oracle formula vs generator truth -----------------------------------

gap <- 0
for (s in tbl$scenario) for (n in ns) {
  d <- generate_scenario_data(s, n, set = "binary", seed = s * 1000 + n)
  oi <- get_oracle_info(s, d$bW, "binary")
  risk <- function(w) {
    eval(parse(text = oi$fmla), envir = c(oi$params, list(X = d$dataset, W = w)))
  }
  gap <- max(gap, abs(risk(1) - d$truth$p1), abs(risk(0) - d$truth$p0))
}
check(gap < 1e-12,
      sprintf("oracle formula reproduces the generator's p0 and p1 (max gap %.1e)", gap))

# ---- 2. tau depends on the described modifiers only ---------------------------

set.seed(2026)
N <- 1e5
cv <- list(X3 = rbinom(N, 1, 0.7), X4 = rnorm(N), X5 = rnorm(N))
tau_at <- function(p, bW, X1, X2) {
  truth_at(p, bW, TRUE, X1, X2,
           if (p$needs_X3) cv$X3, if (p$needs_X4) cv$X4, if (p$needs_X5) cv$X5)$tau
}

shift <- 0
for (s in tbl$scenario) {
  p <- tbl[tbl$scenario == s, ]
  bW <- calibrate_bW(p, 100, "prop")
  a <- tau_at(p, bW, rbinom(N, 1, 0.4), rnorm(N))
  b <- tau_at(p, bW, rbinom(N, 1, 0.4), rnorm(N, 0, 3))
  shift <- max(shift, abs(a - b))
}
check(shift < 1e-12,
      sprintf("tau is unchanged by redrawing X1 and X2 (max shift %.1e)", shift))

# ---- 3 & 4. SD(tau), ATE, power and treated-risk bounds, per scenario x n ------

out <- list()
for (s in tbl$scenario) for (n in ns) {
  p <- tbl[tbl$scenario == s, ]
  bW <- calibrate_bW(p, n, "prop")
  p0 <- control_event_rate(p)
  ate <- bW + te_moments(p)$mean
  b <- treated_risk_bounds(p, bW)
  out[[length(out) + 1]] <- data.frame(
    scenario = s, n = n, bW = bW,
    ate = ate, planned = planned_rd(p, n),
    sd_tau = sd(tau_at(p, bW, rep(0, N), rep(0, N))),
    power = power.prop.test(n = n / 2, p1 = p0, p2 = p0 + ate)$power,
    floor = b[1], ceiling = b[2]
  )
}
res <- do.call(rbind, out)

sd_spread <- tapply(res$sd_tau, res$scenario, function(x) diff(range(x)))
check(sd_spread < 1e-12, "SD(tau) is the same at every n")
check(res$sd_tau[res$scenario == 1] < 1e-12, "scenario 1 is an exact null")
check(abs(res$ate - res$planned) <= ROUNDING + 1e-9,
      "the true ATE is the planned RD, to bW's rounding")
check(abs(res$power - TARGET_POWER) < 0.01, "power is 80% to within 1 point")
check(res$floor >= RD_EPS - ROUNDING & res$ceiling <= 1 - RD_EPS + ROUNDING,
      sprintf("every treated risk is inside [%.2f, %.2f] at n = %s",
              RD_EPS, 1 - RD_EPS, paste(ns, collapse = ", ")))

mb <- t(vapply(miss$scenario[-1], function(s) {
  p <- miss[miss$scenario == s, ]
  treated_risk_bounds(p, calibrate_bW(p, MISSING_N, "prop"), p$bU)
}, numeric(2)))
check(mb[, 1] >= RD_EPS - ROUNDING & mb[, 2] <= 1 - RD_EPS + ROUNDING,
      sprintf("... and under missing/binary's MNAR-Y at n = %d (bU = %.2f)",
              MISSING_N, miss$bU[1]))

# ---- 5. RD_SCALE is the largest feasible scale --------------------------------
# With a unit scale (the continuous coefficients, sign-reversed), g's
# deviations from its mean span [lo, hi], lo <= 0 <= hi. At scale k the treated
# risk's floor is p0_lo + delta_n + k * lo and its ceiling
# p0_hi + delta_n + k * hi (-/+ bU under MNAR-Y), so each bound caps k. A
# negative cap means no k satisfies that bound, which fails the check below.

cts <- resolve_set("continuous")
derived <- lapply(tbl$scenario[-1], function(s) {
  unit <- tbl[tbl$scenario == s, ]
  for (nm in c("b3", "b4", "b34", "b45")) unit[[nm]] <- -cts[[nm]][s]
  dev <- te_range(unit) - te_moments(unit)$mean
  caps <- c()
  for (n in ns) {
    d <- planned_rd(unit, n)
    caps[paste0("floor, n = ", n)] <- (unit$p0_lo + d - RD_EPS) / -dev[1]
    caps[paste0("ceiling, n = ", n)] <- (1 - RD_EPS - unit$p0_hi - d) / dev[2]
  }
  if (s %in% miss$scenario) {
    d <- planned_rd(unit, MISSING_N)
    bU <- abs(miss$bU[1])
    caps["MNAR-Y floor, n = 500"] <- (unit$p0_lo + d - bU - RD_EPS) / -dev[1]
    caps["MNAR-Y ceiling, n = 500"] <- (1 - RD_EPS - unit$p0_hi - d - bU) / dev[2]
  }
  list(k_max = min(caps), binds = names(which.min(caps)))
})
k_max <- vapply(derived, `[[`, numeric(1), "k_max")
check(RD_SCALE[-1] <= k_max + 1e-9 & k_max - RD_SCALE[-1] < 1e-3,
      "RD_SCALE is the largest feasible scale, floored to 3 dp")

# ---- summary ----------------------------------------------------------------

cat("\nper scenario (the ATE at each n is the same in every one):\n")
smry <- data.frame(
  scenario = tbl$scenario,
  RD_SCALE = RD_SCALE,
  k_max = c(NA, k_max),
  binds = c(NA, vapply(derived, `[[`, character(1), "binds")),
  sd_tau = res$sd_tau[res$n == ns[1]],
  floor = tapply(res$floor, res$scenario, min),
  ceiling = tapply(res$ceiling, res$scenario, max)
)
num <- vapply(smry, is.double, logical(1))
smry[num] <- lapply(smry[num], round, 3)
options(width = 200)
print(smry, row.names = FALSE)
