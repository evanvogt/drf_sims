##########
# title: binary study - check the HTE scenarios are what they claim to be
##########
# The binary DGM defines each scenario's treatment effect on the LOGIT scale,
# but the estimand (truth$tau) is a RISK DIFFERENCE. The two only agree about
# which covariates are effect modifiers on the logit scale: on the RD scale the
# prognostic covariates X1 and X2 modify the effect too, through the link, in
# every scenario including "No HTE". This script measures how much.
#
# Checks, per scenario x n (bW is calibrated per n, so the answer moves with n):
#   1. the oracle formula string and the generator's truth agree (the two are
#      separate columns of the scenario table and could drift)
#   2. var(tau_RD) split into the share explained by the described modifiers,
#      E[tau | X3, X4, X5], the share explained by the prognostic covariates,
#      E[tau | X1, X2], and the remainder (their interaction via the link)
#   3. risk range, and the power the marginal risks actually give against the
#      75% bW was calibrated for
#
# Writes nothing. Run from binary/:  Rscript bin_verify_hte.R

source(here::here("R", "dgm_scenarios.R"))

tbl <- resolve_set("binary")
ns <- c(100, 250, 500, 1000)

tau_at <- function(p, bW, X1, X2, X3, X4, X5) {
  truth_at(p, bW, TRUE, X1, X2,
           if (p$needs_X3) X3, if (p$needs_X4) X4, if (p$needs_X5) X5)
}

# ---- 1. oracle formula vs generator truth -----------------------------------

oracle_gap <- 0
for (s in tbl$scenario) for (n in ns) {
  d <- generate_scenario_data(s, n, set = "binary", seed = s * 1000 + n)
  oi <- get_oracle_info(s, d$bW, "binary")
  lp <- function(w) {
    eval(parse(text = oi$fmla), envir = c(oi$params, list(X = d$dataset, W = w)))
  }
  oracle_gap <- max(oracle_gap,
                    abs(plogis(lp(1)) - d$truth$p1),
                    abs(plogis(lp(0)) - d$truth$p0))
}
cat("max |plogis(oracle formula) - generator truth|:",
    format(oracle_gap, digits = 3), "\n\n")

# ---- 2 & 3. where the RD-scale heterogeneity comes from ---------------------

set.seed(2026)
N <- 1e5
cv <- list(X1 = rbinom(N, 1, 0.4), X2 = rnorm(N), X3 = rbinom(N, 1, 0.7),
           X4 = rnorm(N), X5 = rnorm(N))

# X1, X2 integrated out exactly: X1 over {0, 1}, X2 by Gauss-Hermite
K <- length(GH_NODES$x)
P <- list(X1 = rep(c(0, 1), each = K), X2 = rep(sqrt(2) * GH_NODES$x, 2),
          w = c(0.6 * GH_NODES$w, 0.4 * GH_NODES$w))
Nm <- 1e4  # modifier draws used for the conditional expectations

out <- list()
for (s in tbl$scenario) for (n in ns) {
  p <- tbl[tbl$scenario == s, ]
  bW <- calibrate_bW(p, n, "prop")
  tr <- tau_at(p, bW, cv$X1, cv$X2, cv$X3, cv$X4, cv$X5)
  v_tot <- var(tr$tau)

  # E[tau | modifiers], one value per modifier draw
  idx <- rep(seq_len(Nm), each = 2 * K)
  t_m <- tau_at(p, bW, rep(P$X1, Nm), rep(P$X2, Nm),
                cv$X3[idx], cv$X4[idx], cv$X5[idx])$tau
  e_m <- colSums(matrix(t_m * rep(P$w, Nm), nrow = 2 * K))
  v_m <- if (any(p$needs_X3, p$needs_X4, p$needs_X5)) var(e_m) else 0

  # E[tau | X1, X2], one value per quadrature node
  jdx <- rep(seq_len(2 * K), each = Nm)
  t_p <- tau_at(p, bW, P$X1[jdx], P$X2[jdx], rep(cv$X3[1:Nm], 2 * K),
                rep(cv$X4[1:Nm], 2 * K), rep(cv$X5[1:Nm], 2 * K))$tau
  e_p <- colMeans(matrix(t_p, nrow = Nm))
  v_p <- sum(P$w * (e_p - sum(P$w * e_p))^2)

  out[[length(out) + 1]] <- data.frame(
    scenario = s, n = n, bW = bW,
    sd_te_logit = sd(qlogis(tr$p1) - qlogis(tr$p0)),
    ate_rd = mean(tr$tau), sd_tau_rd = sqrt(v_tot),
    share_modifiers = v_m / v_tot,
    share_X1X2 = v_p / v_tot,
    share_link_interaction = 1 - (v_m + v_p) / v_tot,
    p_min = min(tr$p0, tr$p1), p_max = max(tr$p0, tr$p1),
    power_marginal = power.prop.test(n = n / 2, p1 = mean(tr$p0),
                                     p2 = mean(tr$p1))$power
  )
}

res <- do.call(rbind, out)
num <- vapply(res, is.double, logical(1))
res[num] <- lapply(res[num], round, 3)
options(width = 200)
print(res, row.names = FALSE)
