##########
# title: checks for the continuous baseline and bW calibration (bug O)
##########
# Run from anywhere in the repo:
#
#   Rscript R/calibration_check.R          # reference = R/dgm_scenarios.R at HEAD
#   Rscript R/calibration_check.R <ref>    # ...or at another git ref
#
# Checks, in order:
#   1. both continuous tables share one baseline (b0, b1, b2)
#   2. te_moments() agrees with a large Monte Carlo draw
#   3. every continuous scenario at every n a study uses: the ATE has
#      CTS_POWER analytically (to within bW's 2 dp rounding), the mean true tau
#      is -delta, and a simulated Welch t-test on generated data gets close to
#      CTS_POWER
#   4. calibrate_bW() and te_moments() consume no RNG
#   5. draw order: for a fixed seed, W, every covariate, U and the error term
#      match the reference version - only Y, the truth and bW move
#   6. the binary sets are untouched: tables, bW and generated datasets
#      identical to the reference
#
# "Reference" is R/dgm_scenarios.R as committed at <ref>, read with `git show`
# and sourced into its own environment beside the working copy. Once the change
# is committed, pass the commit before it.
#
# Ends by printing the continuous bW table - compare it with a run on the
# cluster: R there is 4.3.2 and here 4.5.3, and a delta sitting on a rounding
# boundary could round to a different bW.

suppressPackageStartupMessages(library(here))

ref <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(ref)) ref <- "HEAD"

new <- new.env()
sys.source(here("R", "dgm_scenarios.R"), envir = new)

old_file <- tempfile(fileext = ".R")
status <- system2("git", c("-C", shQuote(here()), "show",
                           paste0(ref, ":R/dgm_scenarios.R")),
                  stdout = old_file)
if (status != 0) stop("git show ", ref, ":R/dgm_scenarios.R failed")
old <- new.env()
sys.source(old_file, envir = old)

pass <- character()
fail <- character()

report <- function(ok, msg) {
  cat(if (ok) "  PASS  " else "  FAIL  ", msg, "\n", sep = "")
  if (ok) pass <<- c(pass, msg) else fail <<- c(fail, msg)
}

CTS_SETS <- c("continuous", "continuous_missing")
BIN_SETS <- c("binary", "binary_ci", "binary_missing")

# every n a continuous study generates at. validation/continuous/ splits
# n = 1000 at interim_prop 0.25-0.75, so its stages run 250 to 750 in 50s
study_ns <- function(set) {
  if (set == "continuous") sort(unique(c(100, 250, 500, 1000, seq(250, 750, 50))))
  else 500
}

scenario_row <- function(env, set, s) {
  tbl <- env$resolve_set(set)
  tbl[tbl$scenario == s, ]
}

# heterogeneity term g on arbitrary covariate draws, bW = 0 as te_moments uses
g_at <- function(p, X3, X4, X5) {
  eval(parse(text = p$te_expr),
       envir = list(bW = 0, n = length(X3), X3 = X3, X4 = X4, X5 = X5,
                    U_term = 0, b3 = p$b3, b4 = p$b4, b5 = p$b5,
                    b34 = p$b34, b45 = p$b45))
}

# =============================================================================
cat("\n=== 1. one baseline across the continuous scenarios ===\n")

for (set in CTS_SETS) {
  tbl <- new$resolve_set(set)
  same <- all(tbl$b0 == 0.4) && all(tbl$b1 == -0.5) && all(tbl$b2 == 1)
  report(same, sprintf("%s: b0 = 0.4, b1 = -0.5, b2 = 1 in all %d scenarios",
                       set, nrow(tbl)))
}

# =============================================================================
cat("\n=== 2. te_moments() against Monte Carlo ===\n")

set.seed(20260926)
N <- 1e7
mc <- list(X1 = rbinom(N, 1, 0.4), X2 = rnorm(N), X3 = rbinom(N, 1, 0.7),
           X4 = rnorm(N), X5 = rnorm(N))
mc_g <- list()

for (set in CTS_SETS) {
  tbl <- new$resolve_set(set)
  for (s in tbl$scenario) {
    p <- tbl[tbl$scenario == s, ]
    g <- g_at(p, mc$X3, mc$X4, mc$X5)
    if (length(g) == 1) g <- rep(g, N)
    mc_g[[paste(set, s)]] <- c(mean = mean(g), var = var(g))
    tm <- new$te_moments(p)
    ok <- abs(tm$mean - mean(g)) < 0.005 && abs(sqrt(tm$var) - sd(g)) < 0.005
    report(ok, sprintf("%s %d: E[g] %.4f vs MC %.4f, sd(g) %.4f vs MC %.4f",
                       set, s, tm$mean, mean(g), sqrt(tm$var), sd(g)))
  }
}

# =============================================================================
cat("\n=== 3. power and ATE at every study n ===\n")
# analytic power uses the MC moments and a variance written out here, not
# te_moments() or calibrate_bW()'s own arithmetic, so it is an independent check.
# bW is rounded to 2 dp, which moves the ATE by up to 0.005 - at n = 1000 that
# is about 0.02 of power, hence the +-0.02 tolerance rather than +-0.01. The
# ATE check below is the tighter one: the rounding alone bounds it at 0.005

n_truth <- 1e6
sim_reps <- 4000

for (set in CTS_SETS) {
  tbl <- new$resolve_set(set)
  for (s in tbl$scenario) {
    p <- tbl[tbl$scenario == s, ]
    m <- mc_g[[paste(set, s)]]
    var0 <- p$b1^2 * p$X1_prob * (1 - p$X1_prob) + p$b2^2 * p$s2^2 + p$s_err^2
    sd_pooled <- sqrt(var0 + m[["var"]] / 2)

    ns <- study_ns(set)
    if (set == "continuous" && s != 3) ns <- intersect(ns, c(100, 250, 500, 1000))
    pow <- ate_gap <- numeric(length(ns))
    for (i in seq_along(ns)) {
      n <- ns[i]
      bW <- new$calibrate_bW(p, n, "t")
      delta <- power.t.test(n = n / 2, sd = sd_pooled, power = new$CTS_POWER)$delta
      pow[i] <- power.t.test(n = n / 2, delta = abs(bW + m[["mean"]]),
                             sd = sd_pooled)$power
      tau <- new$truth_at(p, bW, FALSE, mc$X1[1:n_truth], mc$X2[1:n_truth],
                          if (p$needs_X3) mc$X3[1:n_truth],
                          if (p$needs_X4) mc$X4[1:n_truth],
                          if (p$needs_X5) mc$X5[1:n_truth])$tau
      ate_gap[i] <- mean(tau) + delta
    }
    report(all(abs(pow - new$CTS_POWER) <= 0.02),
           sprintf("%s %d: analytic power %.3f-%.3f over n = %s", set, s,
                   min(pow), max(pow), paste(range(ns), collapse = "-")))
    report(all(abs(ate_gap) < 0.01),
           sprintf("%s %d: mean true tau within %.4f of -delta", set, s,
                   max(abs(ate_gap))))
  }
}

cat(sprintf("\n  simulated Welch t-test, %d generated datasets per cell\n", sim_reps))
for (set in CTS_SETS) {
  tbl <- new$resolve_set(set)
  sim_ns <- if (set == "continuous") c(100, 250, 500, 1000) else 500
  for (s in tbl$scenario) {
    sim_pow <- vapply(sim_ns, function(n) {
      mean(replicate(sim_reps, {
        d <- new$generate_scenario_data(s, n, set, return_truth = FALSE)$dataset
        t.test(d$Y[d$W == 1], d$Y[d$W == 0])$p.value < 0.05
      }))
    }, numeric(1))
    report(all(abs(sim_pow - new$CTS_POWER) <= 0.03),
           sprintf("%s %d: simulated power %s at n = %s", set, s,
                   paste(sprintf("%.3f", sim_pow), collapse = " / "),
                   paste(sim_ns, collapse = " / ")))
  }
}

# =============================================================================
cat("\n=== 4. no RNG consumed by the calibration ===\n")

set.seed(1)
seed_before <- .Random.seed
for (set in c(CTS_SETS, BIN_SETS)) {
  tbl <- new$resolve_set(set)
  for (s in tbl$scenario) {
    p <- tbl[tbl$scenario == s, ]
    invisible(new$te_moments(p))
    for (n in c(100, 250, 500, 1000)) {
      invisible(new$calibrate_bW(p, n, new$calibration_for(set)))
    }
  }
}
report(identical(.Random.seed, seed_before),
       ".Random.seed unchanged by te_moments() and calibrate_bW() over every set")

# =============================================================================
cat("\n=== 5. draw order unchanged: only Y, the truth and bW move ===\n")

mechs_for <- function(set, s) {
  if (set != "continuous_missing") return(list(NULL))
  m <- list(NULL, "MAR", "MNAR", "MNAR-Y")
  if (s == 1) m[-4] else m
}

# err = Y - p0 - W * (tau + U_term); tau excludes U_term under MNAR-Y
err_of <- function(gen, p, mech) {
  u_term <- if (identical(mech, "MNAR-Y")) p$bU * gen$truth$U else 0
  d <- gen$dataset
  d$Y - gen$truth$p0 - d$W * (gen$truth$tau + u_term)
}

for (set in CTS_SETS) {
  tbl <- new$resolve_set(set)
  for (s in tbl$scenario) {
    for (mech in mechs_for(set, s)) {
      seed <- 1000 * s + 7
      o <- old$generate_scenario_data(s, 250, set, mech = mech, seed = seed)
      w <- new$generate_scenario_data(s, 250, set, mech = mech, seed = seed)
      keep <- setdiff(names(w$dataset), "Y")
      same_x <- identical(o$dataset[keep], w$dataset[keep]) &&
        identical(o$truth$U, w$truth$U)
      same_err <- isTRUE(all.equal(err_of(o, scenario_row(old, set, s), mech),
                                   err_of(w, scenario_row(new, set, s), mech),
                                   tolerance = 1e-12))
      report(same_x && same_err,
             sprintf("%s %d%s: W, covariates, U and error term identical",
                     set, s, if (is.null(mech)) "" else paste0(" (", mech, ")")))
    }
  }
}

# =============================================================================
cat("\n=== 6. binary sets untouched ===\n")

for (set in BIN_SETS) {
  report(identical(old$resolve_set(set), new$resolve_set(set)),
         sprintf("%s: scenario table identical", set))
  tbl <- new$resolve_set(set)
  same_bw <- all(vapply(tbl$scenario, function(s) {
    all(vapply(c(100, 250, 500, 1000), function(n) {
      identical(old$calibrate_bW(scenario_row(old, set, s), n, "prop"),
                new$calibrate_bW(scenario_row(new, set, s), n, "prop"))
    }, logical(1)))
  }, logical(1)))
  report(same_bw, sprintf("%s: bW identical for every scenario and n", set))

  mech_list <- if (set == "binary_missing") list(NULL, "MAR", "MNAR", "MNAR-Y") else list(NULL)
  same_gen <- TRUE
  for (s in tbl$scenario) {
    for (mech in mech_list) {
      if (identical(mech, "MNAR-Y") && s == 1) next
      seed <- 1000 * s + 11
      same_gen <- same_gen && identical(
        old$generate_scenario_data(s, 250, set, mech = mech, seed = seed),
        new$generate_scenario_data(s, 250, set, mech = mech, seed = seed)
      )
    }
  }
  report(same_gen, sprintf("%s: generated datasets and truth identical", set))
}

# =============================================================================
cat("\n=== continuous bW (compare with a run on the cluster) ===\n")

for (set in CTS_SETS) {
  tbl <- new$resolve_set(set)
  ns <- if (set == "continuous") c(100, 250, 500, 1000) else 500
  bw <- t(vapply(tbl$scenario, function(s) {
    vapply(ns, function(n) new$calibrate_bW(scenario_row(new, set, s), n, "t"),
           numeric(1))
  }, numeric(length(ns))))
  if (length(ns) == 1) bw <- t(bw)
  dimnames(bw) <- list(paste0("scenario ", tbl$scenario), paste0("n=", ns))
  cat("\n", set, "\n", sep = "")
  print(bw)
}

# =============================================================================
cat("\n=== summary ===\n")
cat(sprintf("  %d passed, %d failed\n", length(pass), length(fail)))
if (length(fail) > 0) {
  cat("\nfailures:\n")
  for (f in fail) cat("  - ", f, "\n", sep = "")
  quit(status = 1)
}
cat("\nall checks passed.\n")
