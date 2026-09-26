##########
# title: shared data-generating mechanisms
##########
# One implementation of what was four forked DGM files:
#   continuous/cts_dgms.R
#   binary/bin_dgms.R
#   confidence_intervals/continuous/cts_ci_dgms.R  (a verbatim copy of the first,
#                                                   as its own header admitted)
#   confidence_intervals/binary/bin_ci_dgms.R
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
#     W, X1, X2, [X3], [X4], [X5], [U], [err], X01, X02, X03, cats
# X3/X4/X5 are drawn only when the scenario needs them, U only for the MNAR
# mechanisms, and err only for continuous outcomes. R/regression_check.R
# fingerprints the generated dataset precisely to catch a change here.

require(dplyr)

# ---- scenario tables --------------------------------------------------------

# treatment-effect expressions, evaluated with bW, n, X3, X4, X5 and U_term in
# scope alongside the scenario's b* parameters
TE_10 <- c(
  "rep(bW, n)",
  "bW + b3 * X3",
  "bW + b4 * X4",
  "bW + b3 * X3 + b4 * X4",
  "bW + b34 * X3 * X4",
  "bW + b3 * X3 + b4 * X4 + b34 * X3 * X4",
  "bW + b45 * X4 * X5",
  "bW + b3 * X3 + b4 * X4 + b45 * X4 * X5",
  "bW + b4 * cos(X4)",
  "bW + b3 * X3 + b4 * exp(-abs(X4))"     # binary differs here only
)

# oracle formulas. Every table here returns a LINEAR PREDICTOR; for binary
# outcomes the model code applies plogis (oracle_link = "logit"). The binary
# missing-data table shares these strings with the continuous one, so
# missing/binary/bin_miss_models.R must pass "logit" too - it passed "identity"
# until bug M. See the header of R/cate_models.R.
ORACLE_10 <- c(
  "b0 + b1*X$X1 + b2*X$X2 + W*bW",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3)",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b4*X$X4)",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3 + b4*X$X4)",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b34*X$X3*X$X4)",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3 + b4*X$X4 + b34*X$X3*X$X4)",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b45*X$X4*X$X5)",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3 + b4*X$X4 + b45*X$X4*X$X5)",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b4*cos(X$X4))",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3 + b4*exp(-abs(X$X4)))"
)

DESC_10 <- c(
  "No HTE",
  "Simple HTE - binary variable",
  "Simple HTE - continuous variable",
  "Two HTE variables",
  "Continuous-binary interaction (X3*X4)",
  "Single effects + interaction (X3 + X4 + X3*X4)",
  "Continuous-continuous interaction (X4*X5)",
  "Single effects + different interaction (X3 + X4 + X4*X5)",
  "Cosine HTE",
  "Exponential HTE"
)

# the missing-data studies use a reduced, RENUMBERED set of six scenarios.
# scenario k here is NOT scenario k above; the correspondence is
#   1 -> 1 (no HTE), 2 -> 2, 3 -> 4, 4 -> 8, 5 -> 9, 6 -> 3
# Scenario 6 was added after the first five had been run, so it goes at the end
# rather than in main-study order: renumbering would change what the existing
# results/missing/*/scenario_<k>/ directories mean.
# U_term carries the unobserved-confounder contribution under MNAR-Y.
TE_MISS <- c(
  "rep(bW, n)",
  "bW + b3 * X3 + U_term",
  "bW + b3 * X3 + b4 * X4 + U_term",
  "bW + b3 * X3 + b4 * X4 + b45 * X4 * X5 + U_term",
  "bW + b4 * cos(X4) + U_term",
  "bW + b4 * X4 + U_term"
)

ORACLE_MISS <- c(
  "b0 + b1*X$X1 + b2*X$X2 + W*bW",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3)",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3 + b4*X$X4)",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3 + b4*X$X4 + b45*X$X4*X$X5)",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b4*cos(X$X4))",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b4*X$X4)"
)

DESC_MISS <- c(
  "No HTE",
  "Simple HTE - binary variable (X3)",
  "Two HTE variables (X3 + X4)",
  "Single effects + interaction (X3 + X4 + X4*X5)",
  "Non-linear HTE (cos(X4))",
  "Simple HTE - continuous variable (X4)"
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
    b3 = c(NA, 2, NA, 0.3, NA, 2, 2, 2, NA, 0.3),
    b4 = c(NA, NA, -1, -1, NA, 0.5, 0.5, 0.5, 1, 0.1),
    b5 = c(NA, NA, NA, NA, NA, NA, -0.5, -0.5, NA, NA),
    b34 = c(NA, NA, NA, NA, 1, -0.5, NA, NA, NA, NA),
    b45 = c(NA, NA, NA, NA, NA, NA, -0.5, -0.5, NA, NA),
    s2 = 1, s4 = 1, s5 = 1, s_err = 0.5,
    needs_X3 = c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
    needs_X4 = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    needs_X5 = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
    te_expr = TE_10, oracle_expr = ORACLE_10
  ),

  binary = scenario_table(
    scenario = 1:10, description = DESC_10,
    X1_prob = 0.4, X3_prob = 0.7,
    b0 = -0.4, b1 = 0.5, b2 = 0.5,
    b3 = c(NA, -0.4, NA, -0.4, NA, 0.2, 0.2, 0.2, 0.2, 0.2),
    b4 = c(NA, NA, 0.2, 0.3, NA, 0.5, 0.5, 0.5, 0.5, -0.1),
    b5 = c(NA, NA, NA, NA, NA, NA, -0.5, -0.5, NA, NA),
    b34 = c(NA, NA, NA, NA, -0.5, -0.5, NA, NA, NA, NA),
    b45 = c(NA, NA, NA, NA, NA, NA, -0.5, -0.5, NA, NA),
    # the binary DGM drew X2/X4/X5 with a literal sd of 1 rather than via these
    # columns; the values are the same, so one code path serves both outcomes
    s2 = 1, s4 = 1, s5 = 1, s_err = NA,
    needs_X3 = c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
    needs_X4 = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    needs_X5 = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
    te_expr = c(TE_10[1:9], "bW + b4 * exp(X4)"),
    oracle_expr = c(ORACLE_10[1:9],
                    "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b4*exp(X$X4))")
  ),

  continuous_missing = scenario_table(
    scenario = 1:6, description = DESC_MISS,
    X1_prob = 0.4, X3_prob = 0.7,
    # the same shared baseline as the main continuous table (bug O)
    b0 = 0.4, b1 = -0.5, b2 = 1,
    b3 = c(NA, 2, 0.3, 2, NA, NA),
    b4 = c(NA, NA, -1, 0.5, 1, -1),
    b5 = c(NA, NA, NA, -0.5, NA, NA),
    b34 = NA,
    b45 = c(NA, NA, NA, -0.5, NA, NA),
    s2 = 1, s4 = 1, s5 = 1, s_err = 0.5,
    bU = 1, sU = 1,
    needs_X3 = c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE),
    needs_X4 = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE),
    needs_X5 = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
    te_expr = TE_MISS, oracle_expr = ORACLE_MISS
  )

)

# The corrected binary missing-data coefficients. b0/b1/b2 come straight from the
# binary table. b3/b4/b5/b45 are taken from the binary scenario each reduced
# scenario corresponds to (1->1, 2->2, 3->4, 4->8, 5->9) - an inference from the
# scenario descriptions, not something the original code recorded, so worth a
# sanity check before the re-run. Scenario 6 was added later specifically as
# binary scenario 3, so its values are copied from there, not inferred.
SCENARIO_SETS$binary_missing_fixed <- transform(
  SCENARIO_SETS$continuous_missing,
  b0 = -0.4, b1 = 0.5, b2 = 0.5,
  b3 = c(NA, -0.4, -0.4, 0.2, 0.2, NA),
  b4 = c(NA, NA, 0.3, 0.5, 0.5, 0.2),
  b5 = c(NA, NA, NA, -0.5, NA, NA),
  b45 = c(NA, NA, NA, -0.5, NA, NA)
)

# which sets produce a binary outcome
BINARY_SETS <- c("binary", "binary_missing")

#' Resolve a scenario-set name to its table
#'
#' @param set one of names(SCENARIO_SETS), or "binary_ci" for the binary CI study,
#'   or "binary_missing" for the missing/binary study (resolves to the
#'   corrected `binary_missing_fixed` table)
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

#' Mean and variance of a scenario's heterogeneity term g
#'
#' g is the treatment effect with bW = 0 and U_term = 0, so te = bW + g.
#' Evaluated from params$te_expr itself, so it cannot drift from the generator:
#' exactly over X3's two points, by Gauss-Hermite over X4 and X5. Deterministic -
#' it consumes no RNG, so it is safe inside calibrate_bW() (see DRAW ORDER in
#' the file header). calibrate_bW() uses only the mean; the variance is what
#' R/calibration_report.R uses to report each scenario's realised power.
#'
#' @param params one-row scenario params
#' @return list(mean = E[g], var = Var(g))
te_moments <- function(params) {
  axes <- list()
  weights <- list()
  if (params$needs_X3) {
    axes$X3 <- c(0, 1)
    weights$X3 <- c(1 - params$X3_prob, params$X3_prob)
  }
  for (v in c("X4", "X5")) {
    if (params[[paste0("needs_", v)]]) {
      axes[[v]] <- sqrt(2) * params[[sub("X", "s", v)]] * GH_NODES$x
      weights[[v]] <- GH_NODES$w
    }
  }
  if (length(axes) == 0) return(list(mean = 0, var = 0))

  grid <- expand.grid(axes, KEEP.OUT.ATTRS = FALSE)
  w <- Reduce(`*`, expand.grid(weights, KEEP.OUT.ATTRS = FALSE))
  g <- eval(
    parse(text = params$te_expr),
    envir = list(bW = 0, n = nrow(grid), X3 = grid$X3, X4 = grid$X4, X5 = grid$X5,
                 U_term = 0, b3 = params$b3, b4 = params$b4, b5 = params$b5,
                 b34 = params$b34, b45 = params$b45)
  )
  m <- sum(w * g)
  list(mean = m, var = sum(w * (g - m)^2))
}

# ---- generation -------------------------------------------------------------

# power the continuous trials are planned for; the binary branch stays at 0.75
CTS_POWER <- 0.80

#' Calibrate the treatment effect to a fixed power
#'
#' Continuous outcomes: each simulated RCT is planned the way a trial usually
#' is - to detect an ATE, assuming the effect is homogeneous. The planned
#' effect delta gives CTS_POWER in an unadjusted two-sample t-test with n / 2
#' per arm, using the outcome SD with no heterogeneity:
#' sqrt(b1^2 p(1 - p) + b2^2 s2^2 + s_err^2). bW is then set so the true ATE,
#' bW + E[g], equals -delta: the plan gets the average effect right but knows
#' nothing of the heterogeneity around it. Var(g) is left out on purpose, so
#' realised power falls below CTS_POWER as heterogeneity grows (61-80% across
#' scenarios 1-10), just as a real trial planned under homogeneity would. MNAR-Y's
#' U_term is left out for the same reason, which also keeps one bW and one
#' truth per scenario across every missingness mechanism.
#'
#' Before bug O this used sd = s_err + s2 (adding SDs, and ignoring b1 and b2)
#' and set bW rather than the ATE to the planned effect, so the true ATE drifted
#' by E[g] - to the opposite sign in scenarios 2, 6 and 8 - and realised power
#' ran from 3% to 100%.
#'
#' Binary outcomes: a two-proportion test at 75% power on bW at the baseline
#' risk plogis(b0). Unchanged, so still subject to the ATE drift above - see
#' binary/README.md.
#'
#' Neither branch consumes RNG.
calibrate_bW <- function(params, n, calibration = c("t", "prop")) {
  if (match.arg(calibration) == "prop") {
    p1_base <- plogis(params$b0)
    p2 <- power.prop.test(n / 2, p2 = p1_base, power = 0.75)$p1
    round(qlogis(p2) - params$b0, digits = 2)
  } else {
    sd_planned <- sqrt(params$b1^2 * params$X1_prob * (1 - params$X1_prob) +
                         params$b2^2 * params$s2^2 + params$s_err^2)
    delta <- power.t.test(n = n / 2, delta = NULL, sd = sd_planned,
                          power = CTS_POWER)$delta
    round(-delta - te_moments(params)$mean, digits = 2)
  }
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

  X3 <- if (params$needs_X3) rbinom(n, 1, params$X3_prob) else NULL
  X4 <- if (params$needs_X4) rnorm(n, 0, params$s4) else NULL
  X5 <- if (params$needs_X5) rnorm(n, 0, params$s5) else NULL
  U <- if (needs_U) rnorm(n, 0, params$sU) else NULL

  err <- if (!binary) rnorm(n, 0, params$s_err) else NULL

  # the unobserved confounder enters the treatment effect only under MNAR-Y
  U_term <- if (!is.null(mech) && mech == "MNAR-Y") params$bU * U else 0

  treatment_effect <- eval(
    parse(text = params$te_expr),
    envir = list(bW = bW, n = n, X3 = X3, X4 = X4, X5 = X5, U_term = U_term,
                 b3 = params$b3, b4 = params$b4, b5 = params$b5,
                 b34 = params$b34, b45 = params$b45)
  )

  lp <- params$b0 + params$b1 * X1 + params$b2 * X2 + W * treatment_effect
  Y <- if (binary) rbinom(n, 1, plogis(lp)) else lp + err

  # unrelated covariates, always drawn so the fold structure is comparable
  X01 <- rnorm(n, 0, 1)
  X02 <- rnorm(n, 0, 1)
  X03 <- rnorm(n, 0, 1)
  cats <- sample(c("A", "B", "C"), size = n, replace = TRUE, prob = c(0.45, 0.3, 0.25))
  X04 <- as.integer(cats == "A")
  X05 <- as.integer(cats == "B")

  dataset_vars <- list(Y = Y, W = W, X1 = X1, X2 = X2)
  if (params$needs_X3) dataset_vars$X3 <- X3
  if (params$needs_X4) dataset_vars$X4 <- X4
  if (params$needs_X5) dataset_vars$X5 <- X5
  dataset_vars <- c(dataset_vars,
                    list(X01 = X01, X02 = X02, X03 = X03, X04 = X04, X05 = X05))

  result <- list(dataset = as.data.frame(dataset_vars), bW = bW)

  if (return_truth) {
    # the missing-data studies remove U so that tau is the CATE given the
    # observed covariates, averaged over U (U is independent of X). With an
    # identity link that is the U = 0 value, since E[U] = 0; with a logit link
    # it is not, so binary MNAR-Y goes through mnar_y_truth() (bug N)
    reduced <- !is.null(mech)
    link_truth <- binary

    if (!reduced) {
      # the non-MNAR path is shared with build_query_grid_truth() below, via
      # truth_at() - kept as one implementation so the query-grid truth cannot
      # drift from the observed-sample truth
      truth <- truth_at(params, bW, link_truth, X1, X2, X3, X4, X5)
    } else {
      base <- params$b0 + params$b1 * X1 + params$b2 * X2
      if (link_truth && mech == "MNAR-Y") {
        truth <- mnar_y_truth(plogis(base), base + treatment_effect - U_term, params)
      } else {
        if (link_truth) {
          p0 <- plogis(base)
          p1 <- plogis(base + treatment_effect - U_term)
        } else {
          p0 <- base
          p1 <- base + treatment_effect - U_term
        }
        truth <- data.frame(p0 = p0, p1 = p1, tau = p1 - p0)
      }
    }

    if (needs_U) truth$U <- U
    result$truth <- truth
  }

  result
}

#' True p0/p1/tau at an arbitrary set of covariate rows
#'
#' Factored out of generate_scenario_data()'s non-MNAR (reduced == FALSE)
#' truth block so that build_query_grid_truth() below cannot silently diverge
#' from what generate_scenario_data() itself reports as ground truth. The MNAR
#' branch (mech != NULL, subtracting U_term) is NOT reproduced here - it stays
#' inline in generate_scenario_data(), since the query grid is only used by
#' the non-missing CI studies, which never pass mech.
#'
#' @param params one-row scenario params, already subset to `scenario` (as
#'   resolve_set(set) then filtered - see get_oracle_info for the pattern)
#' @param bW calibrated treatment coefficient
#' @param link_truth TRUE to report p0/p1 on the plogis scale (binary outcomes)
#' @param X1,X2 numeric vectors, same length, the two covariates that always
#'   exist
#' @param X3,X4,X5 numeric vectors (same length as X1) or NULL, matching that
#'   scenario's needs_X3/X4/X5 flags
truth_at <- function(params, bW, link_truth, X1, X2, X3 = NULL, X4 = NULL, X5 = NULL) {
  treatment_effect <- eval(
    parse(text = params$te_expr),
    envir = list(bW = bW, n = length(X1), X3 = X3, X4 = X4, X5 = X5, U_term = 0,
                 b3 = params$b3, b4 = params$b4, b5 = params$b5,
                 b34 = params$b34, b45 = params$b45)
  )

  base <- params$b0 + params$b1 * X1 + params$b2 * X2

  if (link_truth) {
    p0 <- plogis(base)
    p1 <- plogis(base + treatment_effect)
  } else {
    p0 <- base
    p1 <- base + treatment_effect
  }

  data.frame(p0 = p0, p1 = p1, tau = p1 - p0)
}

# ---- MNAR-Y truth on the logit scale (bug N) ---------------------------------
# Under MNAR-Y the unobserved U enters the treated arm's linear predictor, so the
# CATE given the observed covariates - what every estimator targets, U being
# unobserved - averages the treated-arm risk over U:
#   p1 = E_U[plogis(eta1 + bU * U)],  U ~ N(0, sU^2),  eta1 = base + te without U
# With an identity link that is eta1 itself, since E[U] = 0, which is why the
# continuous studies need nothing here. With a logit link it is not: the U = 0
# value plogis(eta1) sits further from 0.5 than the average does. The quadrature
# uses GH_NODES, defined above the generation section.

#' E[plogis(eta + s * Z)] for Z ~ N(0, 1), by Gauss-Hermite quadrature
#'
#' Deterministic - it consumes no RNG, so it is safe inside
#' generate_scenario_data() (see DRAW ORDER in the file header).
#'
#' @param eta numeric vector of linear predictors
#' @param s standard deviation of the normal term added to each; 0 returns
#'   plogis(eta) exactly
logistic_normal_mean <- function(eta, s) {
  if (s == 0) return(plogis(eta))
  drop(plogis(outer(eta, sqrt(2) * s * GH_NODES$x, "+")) %*% GH_NODES$w)
}

#' True p0/p1/tau under MNAR-Y for a binary outcome
#'
#' One implementation shared by generate_scenario_data() and
#' repair_mnar_y_truth(), so new runs and repaired old ones cannot disagree.
#' `tau_u0` keeps the old, U = 0 definition - for comparison, and as the marker
#' repair_mnar_y_truth() uses to leave an already-averaged truth alone.
#'
#' @param p0 control-arm risk, plogis(base) - U never enters the control arm
#' @param eta1 treated-arm linear predictor with U removed, base + te
#' @param params one-row scenario params (supplies bU and sU)
mnar_y_truth <- function(p0, eta1, params) {
  p1 <- logistic_normal_mean(eta1, abs(params$bU) * params$sU)
  data.frame(p0 = p0, p1 = p1, tau = p1 - p0, tau_u0 = plogis(eta1) - p0)
}

#' Rebuild binary MNAR-Y truths saved at U = 0, in a collected results tibble
#'
#' Runs made before the bug N fix saved p1 = plogis(eta1), the U = 0 value. p0
#' is still right and eta1 = qlogis(p1), so the averaged truth is recoverable
#' exactly, with no re-run. A truth already carrying `tau_u0` came from the
#' fixed generator and is left alone, so this is safe over a collection holding
#' runs from either side of the fix, and applying it twice changes nothing.
#'
#' @param all_results_df output of get_results(): one row per parameter
#'   combination, `results` a list of list(run, result)
#' @param set the scenario set the study generated from. Supplies bU / sU, and
#'   makes this a no-op for a continuous set.
repair_mnar_y_truth <- function(all_results_df, set) {
  if (!is_binary_set(set) || !"mechanism" %in% names(all_results_df)) {
    return(all_results_df)
  }
  tbl <- resolve_set(set)

  for (i in which(all_results_df$mechanism == "MNAR-Y")) {
    params <- tbl[tbl$scenario == all_results_df$scenario[i], ]
    all_results_df$results[[i]] <- lapply(all_results_df$results[[i]], function(r) {
      truth <- r$result$truth
      if (is.null(truth) || "tau_u0" %in% names(truth)) return(r)
      fixed <- mnar_y_truth(truth$p0, qlogis(truth$p1), params)
      if (!is.null(truth$U)) fixed$U <- truth$U
      # verbatim, not row.names<-, which would turn a row subset's integer
      # row names (complete_cases / IPW) into character ones
      attr(fixed, "row.names") <- attr(truth, "row.names")
      r$result$truth <- fixed
      r
    })
  }
  all_results_df
}

#' Oracle formula and parameter values for a scenario
#'
#' @return list(fmla = <string>, params = <named list>) as run_dr_oracle expects
get_oracle_info <- function(scenario, bW, set) {
  params <- resolve_set(set)
  params <- params[params$scenario == scenario, ]

  param_list <- list(b0 = params$b0, b1 = params$b1, b2 = params$b2, bW = bW)
  for (nm in c("b3", "b4", "b5", "b34", "b45")) {
    v <- params[[nm]]
    if (!is.null(v) && !is.na(v)) param_list[[nm]] <- v
  }

  list(fmla = params$oracle_expr, params = param_list)
}

# ---- covariate query grid (confidence_intervals/{binary,continuous} only) --

# Fixed reference value for covariates the query grid holds constant (X1, X2,
# any of X3/X4/X5 this scenario doesn't need, and the unrelated X01..X05).
# Not a neutral choice for binary scenarios: true tau at a grid point is
# plogis(base + treatment_effect) - plogis(base) where
# base = b0 + b1*X1 + b2*X2, so this reference shifts every binary grid
# point's true CATE (through the nonlinear link), even though it has no
# bearing on how honest the estimators are about hitting whatever that target
# is. For continuous scenarios truth is exactly treatment_effect, which never
# involves X1/X2 at all, so the choice is provably inert there.
GRID_REFERENCE_VALUE <- 0

#' Fixed covariate-grid query points for a scenario's active HTE covariates
#'
#' Varies only the covariates that scenario's treatment effect actually
#' depends on (X3 if needs_X3, X4 if needs_X4, X5 if needs_X5), at fixed
#' design points rather than data-adaptive ones, since every scenario's
#' covariate distributions (X1_prob, X3_prob, s2, s4, s5) are constants, not
#' drawn per replicate - so the same grid is valid, and comparable, across
#' every run of a given scenario. Everything else - X1, X2, any of X3/X4/X5
#' this scenario does not need, and the unrelated X01..X05 - is held at
#' GRID_REFERENCE_VALUE.
#'
#' Scenario 1 ("No HTE") needs none of X3/X4/X5, so the grid degenerates to a
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
  if (isTRUE(params$needs_X3)) active$X3 <- c(0, 1)
  if (isTRUE(params$needs_X4)) active$X4 <- seq(-2, 2, length.out = 5)
  if (isTRUE(params$needs_X5)) active$X5 <- seq(-2, 2, length.out = 5)

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
#'   whichever of X3/X4/X5 that scenario needs
#' @return data.frame(p0, p1, tau), one row per row of grid_df
build_query_grid_truth <- function(scenario, set, bW, grid_df) {
  params <- resolve_set(set)
  params <- params[params$scenario == scenario, ]
  binary <- is_binary_set(set)
  link_truth <- binary

  truth_at(params, bW, link_truth,
           X1 = grid_df$X1, X2 = grid_df$X2,
           X3 = if (isTRUE(params$needs_X3)) grid_df$X3 else NULL,
           X4 = if (isTRUE(params$needs_X4)) grid_df$X4 else NULL,
           X5 = if (isTRUE(params$needs_X5)) grid_df$X5 else NULL)
}
