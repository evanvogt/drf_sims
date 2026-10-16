##########
# title: correlated sample-size sets - SD(tau) and cor(m0, tau) by scenario
##########
# Prints, for both outcomes, both rhos (CORR_RHOS) and scenarios 1-10:
#   sd_tau      SD of the true CATE, exact: sqrt(te_moments()$var), the latent
#               Gauss-Hermite quadrature calibrate_bW() uses (R/dgm_scenarios.R)
#   sd_tau_mc   the same from one large draw of the covariates - a check on
#               sd_tau, which it should match to about 0.005
#   cor_m0_tau  correlation of the control outcome mean m0(x) with tau(x),
#               from that draw. 0 at rho = 0, where the prognostic X1, X2 and
#               the modifiers X3-X5 are independent; NA in scenario 1, where
#               tau is constant
# None of these depend on n: bW only shifts tau. The true ATE, bW and power at
# each n are R/calibration_report.R's.
#
# Needs no simulation output. Writes nothing. Run from sample_size/correlated/:
#   Rscript corr_truth_summary.R

suppressPackageStartupMessages(library(here))
suppressMessages(source(here("R", "dgm_scenarios.R")))

N_MC <- 1e6

rows <- list()
for (outcome in c("continuous", "binary")) {
  binary <- outcome == "binary"
  for (rho in CORR_RHOS) {
    tbl <- resolve_set(corr_set(outcome, rho))
    for (s in tbl$scenario) {
      p <- tbl[tbl$scenario == s, ]
      set.seed(2026)
      z <- correlated_covariates(N_MC, p)
      m0 <- control_mean(p, z$X1, z$X2, binary)
      tau <- truth_at(p, 0, binary, z$X1, z$X2, z$X3, z$X4, z$X5)$tau
      rows[[length(rows) + 1]] <- data.frame(
        outcome = outcome,
        rho = rho,
        scenario = s,
        sd_tau = sqrt(te_moments(p)$var),
        sd_tau_mc = sd(tau),
        cor_m0_tau = if (sd(tau) > 0) cor(m0, tau) else NA_real_
      )
    }
  }
}
res <- do.call(rbind, rows)

num <- c("sd_tau", "sd_tau_mc", "cor_m0_tau")
res[num] <- lapply(res[num], round, 3)
options(width = 200)
for (o in unique(res$outcome)) {
  cat("\n", o, ":\n", sep = "")
  print(res[res$outcome == o, -1], row.names = FALSE)
}
