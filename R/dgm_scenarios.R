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

# ---- flags pending a decision ----------------------------------------------
# All three reproduce known copy-paste bugs so that the move into R/ can be
# proved inert first. Step 8 flips them and deletes the options.

# Bug A: confidence_intervals/binary/bin_ci_dgms.R carried the CONTINUOUS
# coefficient table (b0 = c(0.4, 0.2, ...), b1 = -0.05, b2 = c(2, 2, ...)) on a
# logit scale. It did calibrate bW and compute truth the binary way, so the
# coefficients are the only thing wrong. TRUE keeps them; FALSE uses the binary
# table and forces confidence_intervals/binary + optimal_sf/bin to re-run.
LEGACY_BIN_CI_PARAMS <- FALSE

# missing/binary/bin_miss_dgms.R was forked from the continuous version and only
# half-converted. THREE things came across unchanged, all fixed together because
# they are one mistake and one re-run:
#   1. the continuous coefficient table (as bug A)
#   2. bW calibrated with power.t.test on s_err + s2 - the two-sample t-test for
#      a continuous outcome - instead of power.prop.test
#   3. truth computed as p0 = b0 + b1*X1 + b2*X2 with no plogis, so truth$tau is
#      a LOG-ODDS difference while every estimator targets a RISK DIFFERENCE
# TRUE reproduces all three; FALSE fixes all three and forces missing/binary to
# re-run.
LEGACY_BIN_MISS <- FALSE

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

# oracle formulas. NOTE the two link conventions in this repo: the non-missing
# tables return a LINEAR PREDICTOR and the model code applies plogis for binary
# outcomes (oracle_link = "logit"); the missing tables bake plogis into the
# string (oracle_link = "identity"). See the header of R/cate_models.R.
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

# the missing-data studies use a reduced, RENUMBERED set of five scenarios.
# scenario k here is NOT scenario k above; the correspondence is
#   1 -> 1 (no HTE), 2 -> 2, 3 -> 4, 4 -> 8, 5 -> 9
# U_term carries the unobserved-confounder contribution under MNAR-Y.
TE_5 <- c(
  "rep(bW, n)",
  "bW + b3 * X3 + U_term",
  "bW + b3 * X3 + b4 * X4 + U_term",
  "bW + b3 * X3 + b4 * X4 + b45 * X4 * X5 + U_term",
  "bW + b4 * cos(X4) + U_term"
)

ORACLE_5 <- c(
  "b0 + b1*X$X1 + b2*X$X2 + W*bW",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3)",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3 + b4*X$X4)",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b3*X$X3 + b4*X$X4 + b45*X$X4*X$X5)",
  "b0 + b1*X$X1 + b2*X$X2 + W*(bW + b4*cos(X$X4))"
)

DESC_5 <- c(
  "No HTE",
  "Simple HTE - binary variable (X3)",
  "Two HTE variables (X3 + X4)",
  "Single effects + interaction (X3 + X4 + X4*X5)",
  "Non-linear HTE (cos(X4))"
)

scenario_table <- function(...) data.frame(..., stringsAsFactors = FALSE)

SCENARIO_SETS <- list(

  continuous = scenario_table(
    scenario = 1:10, description = DESC_10,
    X1_prob = 0.4, X3_prob = 0.7,
    b0 = c(0.4, 0.2, 0.3, 0.4, 0.4, 1, 1, 1, 0.4, 0.4),
    b1 = -0.05,
    b2 = c(2, 2, 2, 2, 2, 2, 2, 2, 1, 2),
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
    scenario = 1:5, description = DESC_5,
    X1_prob = 0.4, X3_prob = 0.7,
    b0 = c(0.4, 0.2, 0.4, 1, 0.4),
    b1 = -0.05,
    b2 = c(2, 2, 2, 2, 1),
    b3 = c(NA, 2, 0.3, 2, NA),
    b4 = c(NA, NA, -1, 0.5, 1),
    b5 = c(NA, NA, NA, -0.5, NA),
    b34 = NA,
    b45 = c(NA, NA, NA, -0.5, NA),
    s2 = 1, s4 = 1, s5 = 1, s_err = 0.5,
    bU = 1, sU = 1,
    needs_X3 = c(FALSE, TRUE, TRUE, TRUE, FALSE),
    needs_X4 = c(FALSE, FALSE, TRUE, TRUE, TRUE),
    needs_X5 = c(FALSE, FALSE, FALSE, TRUE, FALSE),
    te_expr = TE_5, oracle_expr = ORACLE_5
  ),

  # LEGACY: identical to continuous_missing, which is the bug. See
  # LEGACY_BIN_MISS_PARAMS and binary_missing_fixed below.
  binary_missing = scenario_table(
    scenario = 1:5, description = DESC_5,
    X1_prob = 0.4, X3_prob = 0.7,
    b0 = c(0.4, 0.2, 0.4, 1, 0.4),
    b1 = -0.05,
    b2 = c(2, 2, 2, 2, 1),
    b3 = c(NA, 2, 0.3, 2, NA),
    b4 = c(NA, NA, -1, 0.5, 1),
    b5 = c(NA, NA, NA, -0.5, NA),
    b34 = NA,
    b45 = c(NA, NA, NA, -0.5, NA),
    # s_err is carried but never used to draw an error term: this study has a
    # binary outcome. It is here because the legacy bW calibration feeds it to
    # power.t.test (see LEGACY_BIN_MISS).
    s2 = 1, s4 = 1, s5 = 1, s_err = 0.5,
    bU = 1, sU = 1,
    needs_X3 = c(FALSE, TRUE, TRUE, TRUE, FALSE),
    needs_X4 = c(FALSE, FALSE, TRUE, TRUE, TRUE),
    needs_X5 = c(FALSE, FALSE, FALSE, TRUE, FALSE),
    te_expr = TE_5,
    oracle_expr = paste0("plogis(", ORACLE_5, ")")
  )
)

# The corrected binary missing-data coefficients. b0/b1/b2 come straight from the
# binary table. b3/b4/b5/b45 are taken from the binary scenario each reduced
# scenario corresponds to (1->1, 2->2, 3->4, 4->8, 5->9) - an inference from the
# scenario descriptions, not something the original code recorded, so worth a
# sanity check before the re-run.
SCENARIO_SETS$binary_missing_fixed <- transform(
  SCENARIO_SETS$binary_missing,
  b0 = -0.4, b1 = 0.5, b2 = 0.5,
  b3 = c(NA, -0.4, -0.4, 0.2, 0.2),
  b4 = c(NA, NA, 0.3, 0.5, 0.5),
  b5 = c(NA, NA, NA, -0.5, NA),
  b45 = c(NA, NA, NA, -0.5, NA)
)

# which sets produce a binary outcome
BINARY_SETS <- c("binary", "binary_missing", "binary_missing_fixed")

#' Resolve a scenario-set name to its table, applying the legacy-bug flags
#'
#' @param set one of names(SCENARIO_SETS), or "binary_ci" for the binary CI study
resolve_set <- function(set) {
  if (set == "binary_ci") {
    # bug A: the CI binary study ran on the continuous table
    return(if (LEGACY_BIN_CI_PARAMS) SCENARIO_SETS$continuous else SCENARIO_SETS$binary)
  }
  if (set == "binary_missing" && !LEGACY_BIN_MISS) {
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
#' Normally this follows the outcome type, but missing/binary inherited the
#' continuous two-sample t-test from the file it was forked from (see
#' LEGACY_BIN_MISS), so it is not derivable from the outcome alone.
calibration_for <- function(set) {
  if (set == "binary_missing" && LEGACY_BIN_MISS) return("t")
  if (is_binary_set(set)) "prop" else "t"
}

# ---- generation -------------------------------------------------------------

#' Calibrate the treatment effect to a fixed power
#'
#' Continuous outcomes use a two-sample t-test, binary outcomes a two-proportion
#' test, both at 75% power. Neither consumes RNG.
calibrate_bW <- function(params, n, calibration = c("t", "prop")) {
  if (match.arg(calibration) == "prop") {
    p1_base <- plogis(params$b0)
    p2 <- power.prop.test(n / 2, p2 = p1_base, power = 0.75)$p1
    round(qlogis(p2) - params$b0, digits = 2)
  } else {
    s_total <- params$s_err + params$s2
    diff <- power.t.test(n = n / 2, delta = NULL, sd = s_total, power = 0.75)$delta
    round(-diff, digits = 2)
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
    base <- params$b0 + params$b1 * X1 + params$b2 * X2

    # the missing-data studies subtract the U contribution so that tau is the
    # marginal treatment effect (U is independent of X, so E[U] = 0)
    reduced <- !is.null(mech)
    # missing/binary reported truth on the linear-predictor scale; see LEGACY_BIN_MISS
    link_truth <- binary && !(set == "binary_missing" && LEGACY_BIN_MISS)

    if (link_truth) {
      p0 <- plogis(base)
      p1 <- plogis(base + treatment_effect - if (reduced) U_term else 0)
    } else {
      p0 <- base
      p1 <- base + treatment_effect - if (reduced) U_term else 0
    }

    truth <- data.frame(p0 = p0, p1 = p1, tau = p1 - p0)
    if (needs_U) truth$U <- U
    result$truth <- truth
  }

  result
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
