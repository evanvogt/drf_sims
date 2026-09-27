##########
# title: shared data-generating mechanisms
##########
# One implementation of what was four forked DGM files:
#   sample_size/continuous/cts_dgms.R
#   sample_size/binary/bin_dgms.R
#   sample_size/confidence_intervals/continuous/cts_ci_dgms.R  (a verbatim copy of the first,
#                                                   as its own header admitted)
#   sample_size/confidence_intervals/binary/bin_ci_dgms.R
#   missing/{continuous,binary}/*_miss_dgms.R      (scenario half; the missingness
#                                                   machinery is in R/missingness.R)
#
# Each fork carried TWO parallel 10-branch switch() statements - one for the
# treatment effect, one for the oracle formula string - that had to be kept in
# step by hand. Both are now string columns on the scenario table, so a scenario
# is defined in one place and the cts/bin divergences are single table cells.
#
# DRAW ORDER IS PART OF THE CONTRACT. Every study seeds with setup_rng_stream()
# and reproduces runs by index, so the sequence of random draws below must not
# change:
#     W, X1, X2, X3, X4, X5, [U], [err], X01, X02, X03, cats
# U only for the MNAR mechanisms, and err only for continuous outcomes.
# R/regression_check.R fingerprints the generated dataset precisely to catch a
# change here.
#
# EVERY SCENARIO DRAWS X1-X5. Since 2026-09-27 X3, X4 and X5 are drawn and
# returned whether or not the scenario's treatment effect uses them, so every
# dataset has the same ten covariates and scenarios differ only in the CATE.
# Before then each was drawn only where the scenario needed it (plus X3 in
# scenario 9), which shifted every later draw: runs of scenarios 1, 2, 4-8 and
# 10 made before 2026-09-27 are not reproducible from this file. Which
# covariates modify the effect is read off te_expr (te_uses()).
#
# SCENARIO NUMBERING. The four scenarios the chapters report come first:
# 1 null, 2 simple, 3 complex, 4 non-linear, then the other six. Before
# 2026-09-26 the ten were numbered differently; old -> new is
#   1 -> 1, 3 -> 2, 8 -> 3, 9 -> 4, 2 -> 5, 4 -> 6, 5 -> 7, 6 -> 8, 7 -> 9, 10 -> 10
# and the missing-data set moved with them (see TE_MISS). Only the numbering
# changed: each scenario generates exactly the data it did under its old number.
# Anything saved before then - collected metrics, figures, results_processing/
# notebooks - uses the old numbers.
#
# BINARY OUTCOMES ARE ON THE RISK-DIFFERENCE SCALE. Every binary set generates
#     P(Y = 1 | x, W) = m0(x) + W * tau(x)
# so the treatment effect adds to the risk, with the control risk m0 bounded in
# [p0_lo, p0_hi] to keep the treated risk inside [0, 1] (control_mean(),
# RD_SCALE). The CATE the binary studies estimate is a risk difference, and it
# is now exactly the scenario's te_expr: X1 and X2 are purely prognostic. Until
# 2026-09-26 the effect was on the logit scale, so the link made X1 and X2
# effect modifiers on the RD scale, scenario 1 was not an RD null, and the HTE
# changed with n (sample_size/binary/README.md). Both outcomes are now one model,
# E[Y | x, W] = m0(x) + W * tau(x), differing in m0, the noise, and the scale and
# sign of the modifiers. Binary scenario 10 now draws X3, which it did not
# before.

require(dplyr)

# ---- scenario tables --------------------------------------------------------

# treatment-effect expressions, evaluated with bW, n, X3, X4, X5 and U_term in
# scope alongside the scenario's b* parameters
TE_10 <- c(
  "rep(bW, n)",
  "bW + b4 * X4",
  "bW + b3 * X3 + b4 * X4 + b45 * X4 * X5",
  "bW + b4 * cos(X4)",
  "bW + b3 * X3",
  "bW + b3 * X3 + b4 * X4",
  "bW + b34 * X3 * X4",
  "bW + b3 * X3 + b4 * X4 + b34 * X3 * X4",
  "bW + b45 * X4 * X5",
  "bW + b3 * X3 + b4 * exp(-abs(X4))"
)

# oracle formulas. Every oracle formula in this file returns the outcome MEAN
# E[Y | X, W] - the linear predictor here, the risk in ORACLE_RD - so the model
# code applies no link. Until the risk-difference DGM the binary formulas were
# linear predictors and the binary studies passed oracle_link = "logit" for the
# model code to apply plogis; a mismatch there was bug M. See the header of
# R/cate_models.R.
ORACLE_10 <- c(
  "b0 + b1*X$X1 + b2*X$X2 + W*bW",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b4*X$X4)",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3 + b4*X$X4 + b45*X$X4*X$X5)",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b4*cos(X$X4))",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3)",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3 + b4*X$X4)",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b34*X$X3*X$X4)",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3 + b4*X$X4 + b34*X$X3*X$X4)",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b45*X$X4*X$X5)",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3 + b4*exp(-abs(X$X4)))"
)

# binary treatment effects, on the RISK-DIFFERENCE scale: TE_10's shapes with X4
# and X5 entering through tanh - bounded, monotone and close to linear within
# +-1 SD - so every treated risk stays in [0, 1] without changing the
# covariates' distributions. cos(X4) and exp(-|X4|) are bounded already.
# Scenario 10 takes the continuous form; it was b4 * exp(X4), unbounded.
TE_RD <- c(
  "rep(bW, n)",
  "bW + b4 * tanh(X4)",
  "bW + b3 * X3 + b4 * tanh(X4) + b45 * tanh(X4) * tanh(X5)",
  "bW + b4 * cos(X4)",
  "bW + b3 * X3",
  "bW + b3 * X3 + b4 * tanh(X4)",
  "bW + b34 * X3 * tanh(X4)",
  "bW + b3 * X3 + b4 * tanh(X4) + b34 * X3 * tanh(X4)",
  "bW + b45 * tanh(X4) * tanh(X5)",
  "bW + b3 * X3 + b4 * exp(-abs(X4))"
)

# binary oracle formulas: the bounded control risk m0 plus W times the treatment
# effect - the risk itself (see control_mean())
RD_CONTROL <- "p0_lo + (p0_hi - p0_lo)*plogis(b0 + b1*X$X1 + b2*X$X2)"
ORACLE_RD <- paste0(RD_CONTROL, c(
  " + W*bW",
  " + W*(bW + b4*tanh(X$X4))",
  " + W*(bW + b3*X$X3 + b4*tanh(X$X4) + b45*tanh(X$X4)*tanh(X$X5))",
  " + W*(bW + b4*cos(X$X4))",
  " + W*(bW + b3*X$X3)",
  " + W*(bW + b3*X$X3 + b4*tanh(X$X4))",
  " + W*(bW + b34*X$X3*tanh(X$X4))",
  " + W*(bW + b3*X$X3 + b4*tanh(X$X4) + b34*X$X3*tanh(X$X4))",
  " + W*(bW + b45*tanh(X$X4)*tanh(X$X5))",
  " + W*(bW + b3*X$X3 + b4*exp(-abs(X$X4)))"
))

# The binary modifiers are the continuous ones with the SIGNS REVERSED, scaled
# onto the risk-difference scale: b = -RD_SCALE[k] * (continuous b). Each
# RD_SCALE[k] is the largest value, floored to 3 dp, that keeps every treated
# risk inside [RD_EPS, 1 - RD_EPS]:
#   - at n = 100, where the planned RD (-0.248) is largest: the floor, which
#     binds in every scenario but 4
#   - at n = 1000, where it is smallest: the ceiling
#   - for scenarios 2-6, also under the missing study's MNAR-Y at n = 500,
#     with bU * tanh(U) added: its ceiling binds scenario 4
# Reversing the signs points each asymmetric scenario's larger swing towards
# LESS benefit, where the bounds leave room - with the continuous signs,
# scenarios 3, 4, 5, 8 and 10 would have to be much smaller. The price is that
# the opposite subgroup benefits more than in continuous/, undoing bug P's sign
# harmonisation. sample_size/binary/bin_verify_hte.R re-derives RD_SCALE from the tables
# and fails if it drifts: change p0_lo, p0_hi, b0-b2, bU, RD_EPS or the binary
# studies' n, and it must be recomputed.
RD_SCALE <- c(NA, 0.082, 0.051, 0.204, 0.137, 0.075, 0.082, 0.137, 0.164, 0.596)
RD_EPS <- 0.01

DESC_10 <- c(
  "No HTE",
  "Simple HTE - continuous variable",
  "Single effects + different interaction (X3 + X4 + X4*X5)",
  "Cosine HTE",
  "Simple HTE - binary variable",
  "Two HTE variables",
  "Continuous-binary interaction (X3*X4)",
  "Single effects + interaction (X3 + X4 + X3*X4)",
  "Continuous-continuous interaction (X4*X5)",
  "Exponential HTE"
)

# the missing-data studies use the first six scenarios above: scenario k here
# is scenario k of the main study. The missing grids run 1-5 (ci_example: 1 and
# 3-6), not all six.
# U_term carries the unobserved-confounder contribution under MNAR-Y.
TE_MISS <- c(
  "rep(bW, n)",
  "bW + b4 * X4 + U_term",
  "bW + b3 * X3 + b4 * X4 + b45 * X4 * X5 + U_term",
  "bW + b4 * cos(X4) + U_term",
  "bW + b3 * X3 + U_term",
  "bW + b3 * X3 + b4 * X4 + U_term"
)

ORACLE_MISS <- c(
  "b0 + b1*X$X1 + b2*X$X2 + W*bW",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b4*X$X4)",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3 + b4*X$X4 + b45*X$X4*X$X5)",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b4*cos(X$X4))",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3)",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3 + b4*X$X4)"
)

DESC_MISS <- c(
  "No HTE",
  "Simple HTE - continuous variable (X4)",
  "Single effects + interaction (X3 + X4 + X4*X5)",
  "Non-linear HTE (cos(X4))",
  "Simple HTE - binary variable (X3)",
  "Two HTE variables (X3 + X4)"
)

scenario_table <- function(...) data.frame(..., stringsAsFactors = FALSE)

SCENARIO_SETS <- list(

  continuous = scenario_table(
    scenario = 1:10, description = DESC_10,
    X1_prob = 0.4, X3_prob = 0.7,
    # one baseline for every scenario (bug O). With the identity link these set
    # the outcome's noise, not the CATE; they varied by scenario before, with
    # b1 = -0.05 leaving X1 all but unprognostic
    b0 = 0.4, b1 = -0.5, b2 = 1,
    # NA wherever the scenario's te_expr doesn't use the coefficient
    b3 = c(NA, NA, 2, NA, 2, 0.3, NA, 2, NA, 0.3),
    b4 = c(NA, -1, 0.5, 1, NA, -1, NA, 0.5, NA, 0.1),
    b34 = c(NA, NA, NA, NA, NA, NA, 1, -0.5, NA, NA),
    b45 = c(NA, NA, -0.5, NA, NA, NA, NA, NA, -0.5, NA),
    s2 = 1, s4 = 1, s5 = 1, s_err = 0.5,
    te_expr = TE_10, oracle_expr = ORACLE_10
  ),

  continuous_missing = scenario_table(
    scenario = 1:6, description = DESC_MISS,
    X1_prob = 0.4, X3_prob = 0.7,
    # the same shared baseline as the main continuous table (bug O)
    b0 = 0.4, b1 = -0.5, b2 = 1,
    b3 = c(NA, NA, 2, NA, 2, 0.3),
    b4 = c(NA, -1, 0.5, 1, NA, -1),
    b34 = NA,
    b45 = c(NA, NA, -0.5, NA, NA, NA),
    s2 = 1, s4 = 1, s5 = 1, s_err = 0.5,
    # under MNAR-Y, U_term = eval(u_expr)
    bU = 1, sU = 1, u_expr = "bU * U",
    te_expr = TE_MISS, oracle_expr = ORACLE_MISS
  )

)

# The binary table is the continuous one on the risk-difference scale (see the
# file header and RD_SCALE): the same covariates and descriptions, the modifiers sign-reversed and scaled by RD_SCALE, and the
# control risk m0 = p0_lo + (p0_hi - p0_lo) * plogis(b0 + b1 X1 + b2 X2). b1 and
# b2 are continuous's, and b0 = -1.72 puts the control event rate at
# E[m0] = 0.400. m0 lies in [0.34, 0.70] (5-95%: 0.35-0.50, SD 0.048), so X1
# and X2 are only weakly prognostic: at n = 100 the planned RD is -0.248, and
# every control risk has to clear that plus the HTE's extra benefit. The event
# is harmful; treatment lowers its risk.
SCENARIO_SETS$binary <- transform(
  SCENARIO_SETS$continuous,
  b0 = -1.72,
  b3 = -RD_SCALE * b3, b4 = -RD_SCALE * b4,
  b34 = -RD_SCALE * b34, b45 = -RD_SCALE * b45,
  s_err = NA,
  te_expr = TE_RD, oracle_expr = ORACLE_RD,
  p0_lo = 0.34, p0_hi = 0.70
)

# The binary missing-data table: the binary table's scenarios 1-6, as scenario
# k of continuous_missing is scenario k of continuous. Under MNAR-Y the
# unobserved U enters the treatment effect as bU * tanh(U): bounded, so the
# treated risk stays inside [RD_EPS, 1 - RD_EPS] (RD_SCALE allows for it), and
# mean zero, so the truth given the observed covariates - the U-free
# m0 + tau - is exactly the average over U, as it is for continuous_missing's
# bU * U. bU = 0.08 gives U's contribution an SD of 0.050, the size of
# scenario 2's HTE. (Until the risk-difference DGM U entered a logit, and that
# average needed quadrature - bug N.)
SCENARIO_SETS$binary_missing_fixed <- transform(
  SCENARIO_SETS$binary[1:6, ],
  description = DESC_MISS,
  bU = 0.08, sU = 1, u_expr = "bU * tanh(U)",
  te_expr = c(TE_RD[1], paste(TE_RD[2:6], "+ U_term"))
)

# which sets produce a binary outcome
BINARY_SETS <- c("binary", "binary_missing")

#' Resolve a scenario-set name to its table
#'
#' @param set one of names(SCENARIO_SETS), or "binary_ci" for the binary CI study,
#'   or "binary_missing" for the missing/binary study (resolves to the
#'   `binary_missing_fixed` table)
resolve_set <- function(set) {
  if (set == "binary_ci") {
    return(SCENARIO_SETS$binary)
  }
  if (set == "binary_missing") {
    return(SCENARIO_SETS$binary_missing_fixed)
  }
  tbl <- SCENARIO_SETS[[set]]
  if (is.null(tbl)) stop("unknown scenario set: ", set)
  tbl
}

is_binary_set <- function(set) {
  set %in% BINARY_SETS || (set == "binary_ci")
}

#' Which power calculation calibrates bW for this set
#'
#' Follows the outcome type: binary outcomes use a two-proportion test,
#' continuous outcomes use a two-sample t-test.
calibration_for <- function(set) {
  if (is_binary_set(set)) "prop" else "t"
}

# ---- quadrature -------------------------------------------------------------

# Gauss-Hermite nodes and weights (Golub-Welsch), built once at source time.
# The weights are normalised to sum to 1, so sum(w * f(sqrt(2) * x)) = E[f(Z)].
GH_NODES <- local({
  k <- 80
  i <- seq_len(k - 1)
  J <- matrix(0, k, k)
  J[cbind(i, i + 1)] <- J[cbind(i + 1, i)] <- sqrt(i / 2)
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = e$vectors[1, ]^2)
})

#' Does a scenario's treatment effect use covariate v?
#'
#' Named in te_expr. Every scenario draws X3-X5, but only the ones te_expr uses
#' get a quadrature, range or query-grid axis.
te_uses <- function(params, v) {
  grepl(paste0("\\b", v, "\\b"), params$te_expr)
}

#' Quadrature grid for a scenario's heterogeneity term g
#'
#' g is the treatment effect with bW = 0 and U_term = 0, so te = bW + g.
#' Evaluated from params$te_expr itself, so it cannot drift from the generator:
#' exactly over X3's two points, by Gauss-Hermite over X4 and X5.
#' Deterministic - it consumes no RNG, so it is safe inside calibrate_bW()
#' (see DRAW ORDER in the file header). Gauss-Hermite is inexact at the kink
#' of scenario 10's exp(-|X4|): E[exp(-|X4|)] comes out 0.5191 against 0.5232,
#' which moves E[g] by 4e-4 at continuous's b4 = 0.1 and 2e-4 at binary's -
#' below bW's rounding either way. RD_SCALE is derived against this E[g], the
#' one the generator uses.
#'
#' @param params one-row scenario params
#' @return list(g = values of g, w = their weights, summing to 1)
te_grid <- function(params) {
  axes <- list()
  weights <- list()
  if (te_uses(params, "X3")) {
    axes$X3 <- c(0, 1)
    weights$X3 <- c(1 - params$X3_prob, params$X3_prob)
  }
  for (v in c("X4", "X5")) {
    if (te_uses(params, v)) {
      axes[[v]] <- sqrt(2) * params[[sub("X", "s", v)]] * GH_NODES$x
      weights[[v]] <- GH_NODES$w
    }
  }

  if (length(axes) == 0) {
    grid <- list()
    w <- 1
  } else {
    grid <- expand.grid(axes, KEEP.OUT.ATTRS = FALSE)
    w <- Reduce(`*`, expand.grid(weights, KEEP.OUT.ATTRS = FALSE))
  }
  g <- eval(
    parse(text = params$te_expr),
    envir = list(bW = 0, n = length(w), X3 = grid$X3, X4 = grid$X4, X5 = grid$X5,
                 U_term = 0, b3 = params$b3, b4 = params$b4,
                 b34 = params$b34, b45 = params$b45)
  )
  list(g = g, w = w)
}

#' Mean and variance of a scenario's heterogeneity term g
#'
#' calibrate_bW() uses the mean; the variance is what R/calibration_report.R
#' uses to report each continuous scenario's realised power.
#'
#' @param params one-row scenario params
#' @return list(mean = E[g], var = Var(g))
te_moments <- function(params) {
  tg <- te_grid(params)
  m <- sum(tg$w * tg$g)
  list(mean = m, var = sum(tg$w * (tg$g - m)^2))
}

#' Range of a scenario's heterogeneity term g over the covariate support
#'
#' X3 at both levels, X4 and X5 over +-6 SD in steps of 0.02 SD - wide and fine
#' enough for every te_expr here: tanh is within 1e-5 of +-1 at 6, a grid point
#' falls within 0.002 of pi for cos(X4), and exp(-|X4|) peaks on the grid at 0.
#' Deterministic, like te_grid().
#'
#' @param params one-row scenario params
#' @return c(min, max) of g
te_range <- function(params) {
  pts <- seq(-6, 6, by = 0.02)
  axes <- list()
  if (te_uses(params, "X3")) axes$X3 <- c(0, 1)
  for (v in c("X4", "X5")) {
    if (te_uses(params, v)) axes[[v]] <- params[[sub("X", "s", v)]] * pts
  }
  grid <- if (length(axes) > 0) expand.grid(axes, KEEP.OUT.ATTRS = FALSE) else list()
  g <- eval(
    parse(text = params$te_expr),
    envir = list(bW = 0, n = max(1, NROW(grid)), X3 = grid$X3, X4 = grid$X4,
                 X5 = grid$X5, U_term = 0, b3 = params$b3, b4 = params$b4,
                 b34 = params$b34, b45 = params$b45)
  )
  range(g)
}

#' Quadrature grid over the prognostic covariates X1 and X2
#'
#' Exact over X1's two points, Gauss-Hermite over X2.
#'
#' @param params one-row scenario params
#' @return list(X1, X2 = covariate values, w = their weights, summing to 1)
baseline_grid <- function(params) {
  x2 <- sqrt(2) * params$s2 * GH_NODES$x
  list(X1 = rep(c(0, 1), each = length(x2)), X2 = c(x2, x2),
       w = c((1 - params$X1_prob) * GH_NODES$w, params$X1_prob * GH_NODES$w))
}

#' Control-arm outcome mean m0(x) = E[Y | x, W = 0]
#'
#' Continuous: the linear predictor b0 + b1 X1 + b2 X2. Binary: that predictor
#' through a logistic scaled into [p0_lo, p0_hi] - the control RISK. Bounding it
#' is what lets the treatment effect add on the risk-difference scale: m0 + tau
#' stays inside [0, 1] however X1 and X2 fall, provided tau stays inside
#' [-p0_lo, 1 - p0_hi] (see RD_SCALE and treated_risk_bounds()).
#'
#' @param params one-row scenario params
#' @param X1,X2 numeric vectors, same length
#' @param binary TRUE for a binary scenario set (is_binary_set())
control_mean <- function(params, X1, X2, binary) {
  eta <- params$b0 + params$b1 * X1 + params$b2 * X2
  if (binary) params$p0_lo + (params$p0_hi - params$p0_lo) * plogis(eta) else eta
}

#' Control event rate E[m0] of a binary scenario
#'
#' Exact over X1, Gauss-Hermite over X2 (baseline_grid()). The p0 the binary
#' power calculation plans from.
#'
#' @param params one-row binary scenario params
control_event_rate <- function(params) {
  base <- baseline_grid(params)
  sum(base$w * control_mean(params, base$X1, base$X2, binary = TRUE))
}

#' Worst-case treated risk over the covariate support, for a binary scenario
#'
#' The treated risk is m0 + bW + g, plus bU * tanh(U) under MNAR-Y. m0 depends
#' only on X1 and X2 and g only on the modifiers, so the extremes add: m0 tends
#' to p0_lo and p0_hi as X2 goes to -/+ infinity, g's range is te_range(), and
#' |tanh(U)| < 1. Every generated risk lies inside these bounds.
#'
#' @param params one-row binary scenario params
#' @param bW calibrated treatment coefficient
#' @param bU the MNAR-Y coefficient; 0 outside MNAR-Y
#' @return c(floor, ceiling)
treated_risk_bounds <- function(params, bW, bU = 0) {
  g <- te_range(params)
  c(params$p0_lo + bW + g[1] - abs(bU), params$p0_hi + bW + g[2] + abs(bU))
}

# ---- generation -------------------------------------------------------------

# power every simulated trial is planned for, continuous and binary
TARGET_POWER <- 0.80

#' Calibrate the treatment effect to a fixed power
#'
#' Each simulated RCT is planned the way a trial usually is - to detect an
#' ATE, assuming the effect is homogeneous - and bW is then set so the true
#' ATE equals the planned effect delta: bW + E[g] = delta. The plan gets the
#' average effect right but knows nothing of the heterogeneity g around it; g
#' enters only through E[g], where it has to. MNAR-Y's U_term is left out, which
#' keeps one bW and one truth per scenario across every missingness mechanism.
#'
#' Continuous outcomes: -delta gives TARGET_POWER in an unadjusted two-sample
#' t-test with n / 2 per arm, using the outcome SD with no heterogeneity:
#' sqrt(b1^2 p(1 - p) + b2^2 s2^2 + s_err^2). Var(g) is left out on purpose, so
#' realised power falls below TARGET_POWER as heterogeneity grows (61-80%
#' across scenarios 1-10), just as a real trial planned under homogeneity would.
#'
#' Binary outcomes: the effect is a risk difference (control_mean()). The plan
#' takes the control event rate p0 = E[m0] and the treated rate p1 that gives
#' TARGET_POWER in a two-proportion test with n / 2 per arm, and
#' delta = p1 - p0. Realised power equals planned: each arm's variance is fixed
#' by its marginal risk, so heterogeneity has no variance to add.
#'
#' Before bug O the continuous branch used sd = s_err + s2 (adding SDs, and
#' ignoring b1 and b2) and set bW rather than the ATE to the planned effect, so
#' the true ATE drifted by E[g] - to the opposite sign in scenarios 3, 5 and 8 -
#' and realised power ran from 3% to 100%. Before bug P the binary branch
#' planned at 75% power at the risk plogis(b0), ignoring X1 and X2, and set bW
#' to the planned log-odds ratio, so the true RD drifted the same way and
#' realised power ran from 5% to 99%. Between bug P and the risk-difference
#' DGM the binary effect was on the logit scale, and bW solved
#' E[plogis(b0 + b1 X1 + b2 X2 + bW + g)] = p1 by uniroot.
#'
#' Neither branch consumes RNG. bW is rounded to 2 dp for a continuous outcome
#' and 3 dp for a binary one: on the risk-difference scale 2 dp would move the
#' RD by up to 0.005, and the power by several percentage points at n = 1000.
calibrate_bW <- function(params, n, calibration = c("t", "prop")) {
  if (match.arg(calibration) == "prop") {
    p0 <- control_event_rate(params)
    delta <- power.prop.test(n = n / 2, p2 = p0, power = TARGET_POWER)$p1 - p0
    digits <- 3
  } else {
    sd_planned <- sqrt(params$b1^2 * params$X1_prob * (1 - params$X1_prob) +
                         params$b2^2 * params$s2^2 + params$s_err^2)
    delta <- -power.t.test(n = n / 2, delta = NULL, sd = sd_planned,
                           power = TARGET_POWER)$delta
    digits <- 2
  }
  round(delta - te_moments(params)$mean, digits = digits)
}

#' Generate one simulated dataset
#'
#' @param scenario scenario index within `set`
#' @param n sample size
#' @param set which scenario table - see SCENARIO_SETS, plus "binary_ci"
#' @param return_truth attach the true p0 / p1 / tau
#' @param mech missingness mechanism; non-NULL draws the unobserved U for the
#'   MNAR variants. "AUX"/"AUX-Y" are accepted as synonyms of "MNAR"/"MNAR-Y",
#'   which is what missing/ci_example still calls them.
#' @param seed optional convenience seed; the studies use setup_rng_stream instead
generate_scenario_data <- function(scenario, n, set, return_truth = TRUE,
                                   mech = NULL, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  params <- resolve_set(set)
  binary <- is_binary_set(set)

  if (!scenario %in% params$scenario) {
    stop("scenario must be one of ", paste(params$scenario, collapse = ", "),
         " for set '", set, "'")
  }
  params <- params[params$scenario == scenario, ]

  # normalise the AUX / MNAR spelling divergence
  if (!is.null(mech)) mech <- sub("^AUX", "MNAR", mech)
  needs_U <- !is.null(mech) && mech %in% c("MNAR", "MNAR-Y")
  if (!is.null(mech) && scenario == 1 && mech == "MNAR-Y") {
    stop("MNAR-Y missingness not applicable to no HTE scenario")
  }

  bW <- calibrate_bW(params, n, calibration_for(set))

  # ---- DRAW ORDER: do not reorder, see the file header ----
  W <- rbinom(n, 1, 0.5)
  X1 <- rbinom(n, 1, params$X1_prob)
  X2 <- rnorm(n, 0, params$s2)

  X3 <- rbinom(n, 1, params$X3_prob)
  X4 <- rnorm(n, 0, params$s4)
  X5 <- rnorm(n, 0, params$s5)
  U <- if (needs_U) rnorm(n, 0, params$sU) else NULL

  err <- if (!binary) rnorm(n, 0, params$s_err) else NULL

  # the unobserved confounder enters the treatment effect only under MNAR-Y, as
  # the set's u_expr: bU * U, or bU * tanh(U) for binary
  U_term <- if (!is.null(mech) && mech == "MNAR-Y") {
    eval(parse(text = params$u_expr), envir = list(bU = params$bU, U = U))
  } else {
    0
  }

  treatment_effect <- eval(
    parse(text = params$te_expr),
    envir = list(bW = bW, n = n, X3 = X3, X4 = X4, X5 = X5, U_term = U_term,
                 b3 = params$b3, b4 = params$b4,
                 b34 = params$b34, b45 = params$b45)
  )

  m0 <- control_mean(params, X1, X2, binary)
  mu <- m0 + W * treatment_effect
  if (binary) {
    # rbinom() returns NA, with only a warning, for a risk outside [0, 1].
    # RD_SCALE keeps every risk inside at the studies' n; this catches any
    # other n, or a table changed without re-deriving it
    if (any(mu < 0 | mu > 1)) {
      stop("scenario ", scenario, " of set '", set, "' at n = ", n,
           " generated a risk outside [0, 1]: range ",
           paste(signif(range(mu), 3), collapse = " to "),
           ". See RD_SCALE and treated_risk_bounds().")
    }
    Y <- rbinom(n, 1, mu)
  } else {
    Y <- mu + err
  }

  # unrelated covariates, always drawn so the fold structure is comparable
  X01 <- rnorm(n, 0, 1)
  X02 <- rnorm(n, 0, 1)
  X03 <- rnorm(n, 0, 1)
  cats <- sample(c("A", "B", "C"), size = n, replace = TRUE, prob = c(0.45, 0.3, 0.25))
  X04 <- as.integer(cats == "A")
  X05 <- as.integer(cats == "B")

  dataset <- data.frame(Y = Y, W = W, X1 = X1, X2 = X2, X3 = X3, X4 = X4, X5 = X5,
                        X01 = X01, X02 = X02, X03 = X03, X04 = X04, X05 = X05)

  result <- list(dataset = dataset, bW = bW)

  if (return_truth) {
    if (is.null(mech)) {
      # shared with build_query_grid_truth() below, via truth_at() - kept as one
      # implementation so the query-grid truth cannot drift from the
      # observed-sample truth
      truth <- truth_at(params, bW, binary, X1, X2, X3, X4, X5)
    } else {
      # the missing-data studies remove U, so that tau is the CATE given the
      # observed covariates, averaged over U (U is independent of X). Every
      # U_term has mean zero - bU * U, or bU * tanh(U) for binary - and enters
      # the outcome mean additively, so that average is the U-free mean on
      # either outcome scale
      p0 <- m0
      p1 <- m0 + treatment_effect - U_term
      truth <- data.frame(p0 = p0, p1 = p1, tau = p1 - p0)
    }

    if (needs_U) truth$U <- U
    result$truth <- truth
  }

  result
}

#' True p0/p1/tau at an arbitrary set of covariate rows
#'
#' Factored out of generate_scenario_data()'s non-MNAR truth block so that
#' build_query_grid_truth() below cannot silently diverge from what
#' generate_scenario_data() itself reports as ground truth. The MNAR branch
#' (mech != NULL, subtracting U_term) is NOT reproduced here - it stays inline
#' in generate_scenario_data(), since the query grid is only used by the
#' non-missing CI studies, which never pass mech.
#'
#' @param params one-row scenario params, already subset to `scenario` (as
#'   resolve_set(set) then filtered - see get_oracle_info for the pattern)
#' @param bW calibrated treatment coefficient
#' @param binary TRUE for a binary set: p0 is then the bounded control risk
#'   (control_mean())
#' @param X1,X2 numeric vectors, same length, the prognostic covariates
#' @param X3,X4,X5 numeric vectors (same length as X1); may be NULL where the
#'   scenario's te_expr does not use them (te_uses())
truth_at <- function(params, bW, binary, X1, X2, X3 = NULL, X4 = NULL, X5 = NULL) {
  treatment_effect <- eval(
    parse(text = params$te_expr),
    envir = list(bW = bW, n = length(X1), X3 = X3, X4 = X4, X5 = X5, U_term = 0,
                 b3 = params$b3, b4 = params$b4,
                 b34 = params$b34, b45 = params$b45)
  )

  p0 <- control_mean(params, X1, X2, binary)
  p1 <- p0 + treatment_effect

  data.frame(p0 = p0, p1 = p1, tau = p1 - p0)
}

#' Oracle formula and parameter values for a scenario
#'
#' @return list(fmla = <string>, params = <named list>) as run_dr_oracle expects
get_oracle_info <- function(scenario, bW, set) {
  params <- resolve_set(set)
  params <- params[params$scenario == scenario, ]

  param_list <- list(b0 = params$b0, b1 = params$b1, b2 = params$b2, bW = bW)
  for (nm in c("b3", "b4", "b34", "b45", "p0_lo", "p0_hi")) {
    v <- params[[nm]]
    if (!is.null(v) && !is.na(v)) param_list[[nm]] <- v
  }

  list(fmla = params$oracle_expr, params = param_list)
}

# ---- covariate query grid (confidence_intervals/{binary,continuous} only) --

# Fixed reference value for covariates the query grid holds constant (X1, X2,
# any of X3/X4/X5 this scenario's te_expr doesn't use, and the unrelated
# X01..X05).
# Inert for the true CATE on either outcome: tau at a grid point is exactly the
# scenario's treatment effect, which never involves X1 or X2. It still sets the
# grid's p0/p1, through m0. (Until the risk-difference DGM the binary CATE was
# plogis(base + te) - plogis(base), so this reference shifted every binary grid
# point's true CATE through the link.)
GRID_REFERENCE_VALUE <- 0

#' Fixed covariate-grid query points for a scenario's active HTE covariates
#'
#' Varies only the covariates that scenario's treatment effect actually
#' depends on (whichever of X3/X4/X5 te_uses() finds in te_expr), at fixed
#' design points rather than data-adaptive ones, since every scenario's
#' covariate distributions (X1_prob, X3_prob, s2, s4, s5) are constants, not
#' drawn per replicate - so the same grid is valid, and comparable, across
#' every run of a given scenario. Everything else - X1, X2, any of X3/X4/X5
#' this scenario's effect does not use, and the unrelated X01..X05 - is held at
#' GRID_REFERENCE_VALUE.
#'
#' Scenario 1 ("No HTE") uses none of X3/X4/X5, so the grid degenerates to a
#' single row at the reference point - handled explicitly (expand.grid() over
#' zero variables returns a 1-row, 0-column data.frame, which is not a useful
#' thing to cbind against) rather than left to that edge case.
#'
#' @param scenario scenario index within `set`
#' @param set which scenario table, as get_oracle_info (use "binary_ci" for
#'   the binary CI study, "continuous" for the continuous one)
#' @param covariate_names names(data)[-(1:2)] from the observed dataset - the
#'   exact column set/order predict(forest, newdata=...) must be handed,
#'   since grf matches newdata columns by position, not name
#' @return data.frame with columns exactly `covariate_names`, in that order
build_query_grid <- function(scenario, set, covariate_names) {
  params <- resolve_set(set)
  params <- params[params$scenario == scenario, ]

  active <- list()
  if (te_uses(params, "X3")) active$X3 <- c(0, 1)
  if (te_uses(params, "X4")) active$X4 <- seq(-2, 2, length.out = 5)
  if (te_uses(params, "X5")) active$X5 <- seq(-2, 2, length.out = 5)

  grid <- if (length(active) > 0) {
    do.call(expand.grid, c(active, list(stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)))
  } else {
    data.frame(row.names = 1)
  }

  reference_names <- setdiff(covariate_names, names(grid))
  for (nm in reference_names) grid[[nm]] <- GRID_REFERENCE_VALUE

  grid[, covariate_names, drop = FALSE]
}

#' True p0/p1/tau at each row of a query grid
#'
#' @param scenario,set,bW as get_oracle_info
#' @param grid_df output of build_query_grid() - must carry X1, X2 and
#'   whichever of X3/X4/X5 that scenario's te_expr uses
#' @return data.frame(p0, p1, tau), one row per row of grid_df
build_query_grid_truth <- function(scenario, set, bW, grid_df) {
  params <- resolve_set(set)
  params <- params[params$scenario == scenario, ]

  truth_at(params, bW, is_binary_set(set),
           X1 = grid_df$X1, X2 = grid_df$X2,
           X3 = grid_df$X3, X4 = grid_df$X4, X5 = grid_df$X5)
}
