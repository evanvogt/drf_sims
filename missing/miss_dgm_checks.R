##########
# title: missing-data DGM checks - does each mechanism do what ADEMP.md says?
##########
# Checks the claims missing/ADEMP.md makes about the DGM, at a large n where the
# estimators' noise is gone. For both outcomes (continuous_missing,
# binary_missing), one scenario, each mechanism:
#   1. covariates: prevalences and correlations of the copula
#   2. amputation: proportion incomplete, missing rate per covariate, patterns
#   3. selection: how well incompleteness is predicted from X1-X5, from the
#      always-observed X01-X05, from the IPW model (X01-X05 + W * Y, since
#      2026-09-30), and from Y in each arm; the IPW weights' spread; how far
#      the missing X4 values sit from the observed; E[U_term | complete /
#      incomplete] and U's share of Var(Y)
#   4. pairing: within a run, MAR and MNAR-Y0 share W, X, err and cats, and
#      differ in Y by U's term alone (R/dgm_scenarios.R, DRAW ORDER)
#   5. mechanism strength on a common scale, both outcomes side by side
#   6. each handling method's signature: a correctly specified
#      lm(Y ~ W * covariates) fit on the handled data, its CATE scored on the
#      complete and the incomplete units separately (bias, slope on the truth,
#      RMSE). IPW is that lm on the complete cases, weighted. Multiple
#      imputation and missForest at a smaller n, as they are slow.
#
# Writes nothing. Run from the repo root:
#   Rscript missing/miss_dgm_checks.R [scenario] [n]
# scenario defaults to 2 (simple HTE on X4), n to 20000. Takes a few minutes,
# most of it the imputation block.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
})
suppressMessages(source(here("R", "missingness.R")))

args <- commandArgs(trailingOnly = TRUE)
SCENARIO <- if (length(args) >= 1) as.integer(args[1]) else 2
N <- if (length(args) >= 2) as.integer(args[2]) else 20000
STUDY_N <- 500        # the studies' n: bW, so the true ATE, is calibrated there
N_IMP_CHECK <- 3000   # n for the multiple imputation / missForest signatures
M_IMP_CHECK <- 5      # imputations there (the study uses 50)
SETS <- c("continuous_missing", "binary_missing")
TYPE <- "both"
PROP <- 0.3
options(width = 200)

# Mann-Whitney AUC of `score` for the 0/1 `y`
auc <- function(score, y) {
  r <- rank(score)
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# the CATE of a correctly specified linear model fit to (handled) data: the
# scenario's modifiers enter through the set's own transform (tanh for binary);
# `w` optional observation weights (the IPW arm)
cate_lm <- function(d, binary, w = NULL) {
  d$f4 <- if (binary) tanh(d$X4) else d$X4
  covs <- c("X1", "X2", "X3", "f4", "X5", "X01", "X02", "X03",
            grep("_missing$", names(d), value = TRUE))
  fit <- lm(as.formula(paste("Y ~ W * (", paste(covs, collapse = " + "), ")")), data = d,
            weights = w)
  d1 <- d
  d1$W <- 1
  d0 <- d
  d0$W <- 0
  predict(fit, d1) - predict(fit, d0)
}

# one handling method's signature, split by unit completeness
signature <- function(label, tau_hat, tau, incomplete) {
  one <- function(keep) {
    if (!any(keep)) return(c(bias = NA, slope = NA, rmse = NA))
    e <- tau_hat[keep] - tau[keep]
    slope <- if (sd(tau[keep]) > 0) unname(coef(lm(tau_hat[keep] ~ tau[keep]))[2]) else NA
    c(bias = mean(e), slope = slope, rmse = sqrt(mean(e^2)))
  }
  a <- one(rep(TRUE, length(tau)))
  cu <- one(!incomplete)
  iu <- one(incomplete)
  data.frame(method = label, n = length(tau),
             bias = a[["bias"]], slope = a[["slope"]], rmse = a[["rmse"]],
             bias_cu = cu[["bias"]], rmse_cu = cu[["rmse"]],
             bias_iu = iu[["bias"]], rmse_iu = iu[["rmse"]])
}

strength <- list()

for (set in SETS) {
  binary <- is_binary_set(set)
  params <- resolve_set(set)
  p <- params[params$scenario == SCENARIO, ]

  cat("\n#################", set, "- scenario", SCENARIO, "-", p$description,
      "- n =", N, "#################\n")

  runs <- list()
  for (mech in MISS_MECHS) {
    if (SCENARIO == 1 && mech == "MNAR-tau") next
    set.seed(SCENARIO * 1000 + 1)   # the same seed for every mechanism: paired
    g <- generate_scenario_data(SCENARIO, N, set, mech = mech, calib_n = STUDY_N)
    runs[[mech]] <- g
  }

  # ---- 1. covariates ----
  cd <- runs[["MAR"]]$dataset
  cat("\n1. covariates\n")
  cat(sprintf("   P(X1 = 1) = %.3f (planned %.2f), P(X3 = 1) = %.3f (planned %.2f)\n",
              mean(cd$X1), p$X1_prob, mean(cd$X3), p$X3_prob))
  cr <- cor(cd[, c("X1", "X2", "X3", "X4", "X5", "X01", "X02", "X03")])
  cont <- c("X2", "X4", "X5", "X01", "X02", "X03")
  cat(sprintf("   mean correlation among the continuous copula columns %.3f (latent rho %.2f); involving X1/X3 %.3f\n",
              mean(cr[cont, cont][upper.tri(cr[cont, cont])]), p$rho,
              mean(c(cr["X1", setdiff(rownames(cr), "X1")], cr["X3", setdiff(rownames(cr), c("X1", "X3"))]))))
  cat(sprintf("   X04/X05 vs the copula: max |correlation| %.3f (drawn independently)\n",
              max(abs(cor(cd[, c("X04", "X05")], cd[, c("X1", "X2", "X3", "X4", "X5")])))))

  # ---- 4. pairing across mechanisms ----
  cat("\n4. pairing (same run, different mechanism)\n")
  same_x <- all.equal(runs[["MAR"]]$dataset[, -1], runs[["MNAR-Y0"]]$dataset[, -1])
  u_term <- eval(parse(text = p$u_expr), envir = list(bU = p$bU, U = runs[["MNAR-Y0"]]$truth$U))
  y_gap <- max(abs(runs[["MNAR-Y0"]]$dataset$Y - runs[["MAR"]]$dataset$Y - u_term))
  cat(sprintf("   MAR vs MNAR-Y0: W and every covariate identical: %s\n", isTRUE(same_x)))
  if (!binary) {
    cat(sprintf("   continuous Y differs by exactly U's term: max gap %.1e\n", y_gap))
  } else {
    cat("   binary Y is redrawn from the shifted risk, so it is paired on the risk, not on Y\n")
  }

  for (mech in names(runs)) {
    g <- runs[[mech]]
    cd <- g$dataset
    U <- if (mech %in% MNAR_MECHS) g$truth$U else NULL
    amp <- suppressWarnings(introduce_missingness(cd, TYPE, PROP, mech, U = U))
    mask <- amputation_mask(amp)
    incomplete <- rowSums(mask) > 0
    inc <- as.integer(incomplete)
    tau <- g$truth$tau

    cat("\n=====", set, "-", mech, "=====\n")

    # ---- 2. amputation ----
    cat("2. amputation\n")
    cat(sprintf("   incomplete %.3f (planned %.2f); missing per covariate: %s\n",
                mean(incomplete), PROP,
                paste(sprintf("%s %.3f", colnames(mask), colMeans(mask)), collapse = ", ")))
    patt <- table(apply(mask[incomplete, , drop = FALSE] * 1, 1, paste, collapse = ""))
    cat(sprintf("   %d distinct missingness patterns among the incomplete units (%d-%d units each)\n",
                length(patt), min(patt), max(patt)))

    # ---- 3. selection ----
    cat("3. selection\n")
    auc_x <- auc(fitted(glm(inc ~ X1 + X2 + X3 + X4 + X5, binomial, cd)), inc)
    auc_aux <- auc(fitted(glm(inc ~ X01 + X02 + X03 + X04 + X05, binomial, cd)), inc)
    auc_ipw <- auc(fitted(glm(inc ~ X01 + X02 + X03 + X04 + X05 + W * Y, binomial, cd)), inc)
    auc_y <- auc(cd$Y, inc)
    auc_y1 <- auc(cd$Y[cd$W == 1], inc[cd$W == 1])
    auc_y0 <- auc(cd$Y[cd$W == 0], inc[cd$W == 0])
    cat(sprintf("   AUC of incompleteness: from X1-X5 %.3f, from X01-X05 %.3f, from the IPW model (X01-X05 + W * Y) %.3f; from Y %.3f (treated %.3f, control %.3f)\n",
                auc_x, auc_aux, auc_ipw, auc_y, auc_y1, auc_y0))
    h_ipw <- suppressMessages(handle_missingness(amp, "IPW"))
    cat(sprintf("   IPW weights (stabilised): mean %.3f, CV %.3f, min %.2f, max %.2f\n",
                mean(h_ipw$ipw), sd(h_ipw$ipw) / mean(h_ipw$ipw),
                min(h_ipw$ipw), max(h_ipw$ipw)))
    miss4 <- mask[, "X4"]
    cat(sprintf("   missing X4 values sit %.2f SD above the observed ones\n",
                (mean(cd$X4[miss4]) - mean(cd$X4[!miss4])) / sd(cd$X4)))
    shift <- c(complete = NA_real_, incomplete = NA_real_)
    var_share <- 0
    if (!is.null(U)) {
      ut <- eval(parse(text = p$u_expr), envir = list(bU = p$bU, U = U))
      shift <- c(complete = mean(ut[!incomplete]), incomplete = mean(ut[incomplete]))
      var_share <- var(ut) / var(cd$Y)
      cat(sprintf("   E[U_term | complete] %.3f, E[U_term | incomplete] %.3f; U_term is %.1f%% of Var(Y)\n",
                  shift[["complete"]], shift[["incomplete"]], 100 * var_share))
    }

    # ---- 6. signatures ----
    sig <- list(signature("complete_data", cate_lm(cd, binary), tau, incomplete))
    tau_cc <- rep(NA_real_, N)
    tau_cc[!incomplete] <- cate_lm(cd[!incomplete, ], binary)
    sig[["cc"]] <- signature("complete_cases", tau_cc[!incomplete], tau[!incomplete],
                             incomplete[!incomplete])
    sig[["ipw"]] <- signature("IPW", cate_lm(h_ipw$data, binary, w = h_ipw$ipw),
                              tau[!incomplete], incomplete[!incomplete])
    for (m in c("mean_imputation", "missing_indicator", "regression")) {
      h <- suppressMessages(handle_missingness(amp, m))
      sig[[m]] <- signature(m, cate_lm(h$data, binary), tau, incomplete)
    }
    sub <- seq_len(N_IMP_CHECK)
    amp_s <- amp[sub, ]
    sig[["cd_s"]] <- signature(sprintf("complete_data (n = %d)", N_IMP_CHECK),
                               cate_lm(cd[sub, ], binary), tau[sub], incomplete[sub])
    # capture.output: missForest and mice print every iteration
    invisible(capture.output(
      mf <- suppressMessages(handle_missingness(amp_s, "missforest"))))
    sig[["mf"]] <- signature(sprintf("missforest (n = %d)", N_IMP_CHECK),
                             cate_lm(mf$data, binary), tau[sub], incomplete[sub])
    invisible(capture.output(
      mi <- suppressMessages(handle_missingness(amp_s, "multiple_imputation",
                                                n_imp = M_IMP_CHECK))))
    sig[["mi"]] <- signature(sprintf("multiple_imputation (n = %d, m = %d)", N_IMP_CHECK, M_IMP_CHECK),
                             rowMeans(sapply(mi$data, cate_lm, binary = binary)),
                             tau[sub], incomplete[sub])
    cat("6. signatures: CATE of a correctly specified lm(Y ~ W * covariates) on the handled data\n")
    sig <- do.call(rbind, sig)
    num <- vapply(sig, is.double, logical(1))
    sig[num] <- lapply(sig[num], signif, 3)
    print(sig, row.names = FALSE)

    # ---- 5. strength, collected for the common-scale table ----
    ate <- mean(tau)
    cc_ate_bias <- mean(tau_cc[!incomplete]) - mean(tau[!incomplete])
    strength[[length(strength) + 1]] <- data.frame(
      set = set, mechanism = mech,
      auc_incomplete_X = auc_x, auc_incomplete_Y = auc_y,
      u_share_varY = var_share,
      cc_ate_bias = cc_ate_bias,
      cc_ate_bias_pct_ate = 100 * cc_ate_bias / abs(ate),
      cc_ate_bias_pct_sd_tau = if (sd(tau) > 0) 100 * cc_ate_bias / sd(tau) else NA,
      shift_incomplete = shift[["incomplete"]]
    )
  }
}

cat("\n5. mechanism strength on a common scale (scenario", SCENARIO, ")\n")
cat("   cc_ate_bias: complete cases' ATE bias from the correctly specified lm, which is\n")
cat("   noise at this n except under MNAR-tau\n")
strength <- do.call(rbind, strength)
num <- vapply(strength, is.double, logical(1))
strength[num] <- lapply(strength[num], signif, 3)
print(strength, row.names = FALSE)
