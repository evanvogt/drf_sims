##########
# title: what the competing-risk DGM implies, at each covariate correlation
##########
# Produces the numbers in ADEMP.md "What the DGM implies" and the separation
# table under "Why these values", for every rho in CORR_RHOS:
#
#   1. event mix by the horizon, per arm - generate_surv_data() at n = 20,000
#      per (rho, scenario, censoring) cell, seeded once (illustrative, not the
#      study's streams);
#   2. population truths - surv_truth() on an X1 x X3 x X2 grid, interpolated
#      in X2 over N_POP copula draws per rho: mean and SD of each tau,
#      cor(tau_RMTL1, tau_RMTL2), control-arm levels;
#   3. the X3 separation on the RMTL scale, two ways:
#        marginal      E[tau | X3 = 1] - E[tau | X3 = 0]. Under correlation
#                      this also picks up the X1/X2-driven heterogeneity of
#                      the units who happen to have X3 = 1;
#        standardised  E[tau(X1, X2, 1) - tau(X1, X2, 0)], X3 switched with
#                      (X1, X2) held at their own distribution - X3's own
#                      effect.
#      At rho = 0 the two agree. separation = smallest labelled X3 difference
#      / largest induced one (RMTL1: labelled 3, 4, 7, induced 5, 6; RMTL2:
#      labelled 5, 6, 7, induced 3, 4).
#
# Local, about 5-10 minutes:  Rscript surv_dgm_check.R  (from competing_risk/)

library(here)
library(dplyr)
library(future)
library(knitr)
source(here("R", "utils.R"))
source(here("competing_risk", "surv_dgm.R"))

plan(multisession, workers = min(6, max(1, parallelly::availableCores() - 1)))

N_MIX <- 20000   # units per event-mix cell
N_POP <- 200000  # copula draws per rho for the population truths
X2_GRID <- seq(-4.5, 4.5, length.out = 41)
horizon <- survival_scenario_params$event_horizon[1]
scenarios <- survival_scenario_params$scenario

tau_cols <- c("tau_RMTL1", "tau_RMTL2", "tau_RMSTc", "tau_RMST1_cs", "tau_RMST2_cs")
ctrl_cols <- c("RMTL1_0", "RMTL2_0", "RMSTc_0", "RMST1_cs_0", "RMST2_cs_0")

show <- function(df, digits = 2) {
  print(kable(df, format = "pipe", digits = digits))
  cat("\n")
}

# ---- 1. event mix -------------------------------------------------------------
cat("=== 1. Event mix by the horizon (n =", N_MIX, "per cell) ===\n\n")
set.seed(20261001)
mix <- bind_rows(lapply(CORR_RHOS, function(rho) {
  bind_rows(lapply(scenarios, function(s) {
    bind_rows(lapply(c(FALSE, TRUE), function(cens) {
      d <- generate_surv_data(s, N_MIX, rho = rho, return_truth = FALSE,
                              censoring = cens)$dataset
      d %>%
        group_by(arm = ifelse(W == 1, "treated", "control")) %>%
        summarise(E1 = mean(D == 1 & Y <= horizon),
                  E2 = mean(D == 2 & Y <= horizon),
                  past = mean(Y > horizon),
                  cens_before = mean(D == 0 & Y < horizon),
                  .groups = "drop") %>%
        mutate(rho = rho, scenario = s, censoring = cens, .before = 1)
    }))
  }))
}))

for (r in CORR_RHOS) {
  cat("rho =", r, "- treated arm; E1 / E2 / event-free without censoring,",
      "censored before / observed past the horizon with it\n\n")
  show(mix %>%
         filter(rho == r, arm == "treated") %>%
         group_by(scenario) %>%
         summarise(E1 = E1[!censoring], E2 = E2[!censoring],
                   event_free = past[!censoring],
                   cens_before_28 = cens_before[censoring],
                   past_28_cens = past[censoring], .groups = "drop"), 3)
  cat("rho =", r, "- control arm, range over scenarios\n\n")
  show(mix %>%
         filter(rho == r, arm == "control") %>%
         group_by(censoring) %>%
         summarise(across(c(E1, E2, past, cens_before),
                          ~ sprintf("%.2f-%.2f", min(.x), max(.x))),
                   .groups = "drop"))
}

# ---- 2. population truths -----------------------------------------------------
cat("=== 2. Population truths (", N_POP, " copula draws per rho) ===\n\n", sep = "")

grid <- expand.grid(X2 = X2_GRID, X1 = 0:1, X3 = 0:1)
grid_truth <- lapply(scenarios, function(s) {
  params <- survival_scenario_params[survival_scenario_params$scenario == s, ]
  bind_cols(grid, surv_truth(params, grid$X1, grid$X2, grid$X3))
})

#' A truth column at each draw: interpolated in X2 within its (X1, X3) cell
at_draws <- function(gt, col, X1, X2, X3) {
  out <- numeric(length(X2))
  for (a in 0:1) for (b in 0:1) {
    g <- gt[gt$X1 == a & gt$X3 == b, ]
    idx <- X1 == a & X3 == b
    out[idx] <- approx(g$X2, g[[col]], X2[idx], rule = 2)$y
  }
  out
}

pop <- list()
sep <- list()
for (r in CORR_RHOS) {
  set.seed(1)
  cv <- correlated_covariates(N_POP, list(X1_prob = 0.4, X3_prob = 0.7, rho = r,
                                          s2 = 1, s4 = 1, s5 = 1))
  X1 <- cv$X1; X2 <- cv$X2; X3 <- cv$X3

  for (s in scenarios) {
    gt <- grid_truth[[s]]
    tau <- sapply(tau_cols, function(cl) at_draws(gt, cl, X1, X2, X3))
    pop[[length(pop) + 1]] <- tibble(
      rho = r, scenario = s,
      !!!setNames(lapply(tau_cols, function(cl) {
        sprintf("%.2f (%.2f)", mean(tau[, cl]), sd(tau[, cl]))
      }), tau_cols),
      cor_RMTL1_RMTL2 = if (sd(tau[, "tau_RMTL1"]) > 0 && sd(tau[, "tau_RMTL2"]) > 0)
        round(cor(tau[, "tau_RMTL1"], tau[, "tau_RMTL2"]), 2) else NA,
      RMSTc_X3_0 = mean(tau[X3 == 0, "tau_RMSTc"]),
      RMSTc_X3_1 = mean(tau[X3 == 1, "tau_RMSTc"])
    )
    for (cl in c("tau_RMTL1", "tau_RMTL2")) {
      std <- mean(at_draws(gt, cl, X1, X2, rep(1, N_POP)) -
                    at_draws(gt, cl, X1, X2, rep(0, N_POP)))
      sep[[length(sep) + 1]] <- tibble(
        rho = r, scenario = s, target = cl,
        marginal = mean(tau[X3 == 1, cl]) - mean(tau[X3 == 0, cl]),
        standardised = std
      )
    }
  }

  if (r == CORR_RHOS[1]) {
    gt <- grid_truth[[1]]
    ctrl <- sapply(ctrl_cols, function(cl) mean(at_draws(gt, cl, X1, X2, X3)))
    cat("control-arm levels (rho =", r, "):\n")
    print(round(ctrl, 1))
    cat("\n")
  }
}
pop <- bind_rows(pop)
sep <- bind_rows(sep)

for (r in CORR_RHOS) {
  cat("rho =", r, "- tau, population mean (SD), in days\n\n")
  show(filter(pop, rho == r) %>% select(-rho, -RMSTc_X3_0, -RMSTc_X3_1))
}
cat("scenario 6, tau_RMSTc by X3 (marginal):\n\n")
show(filter(pop, scenario == 6) %>% select(rho, RMSTc_X3_0, RMSTc_X3_1))

# ---- 3. separation ------------------------------------------------------------
cat("=== 3. X3 separation on the RMTL scale ===\n\n")
labelled <- list(tau_RMTL1 = c(3, 4, 7), tau_RMTL2 = c(5, 6, 7))
induced  <- list(tau_RMTL1 = c(5, 6),    tau_RMTL2 = c(3, 4))

show(sep %>%
       filter(scenario %in% 3:7) %>%
       tidyr::pivot_wider(names_from = target,
                          values_from = c(marginal, standardised)) %>%
       arrange(rho, scenario))

separation <- sep %>%
  group_by(rho, target) %>%
  summarise(
    across(c(marginal, standardised), function(v) {
      lab <- abs(v[scenario %in% labelled[[cur_group()$target]]])
      ind <- abs(v[scenario %in% induced[[cur_group()$target]]])
      min(lab) / max(ind)
    }),
    .groups = "drop"
  )
cat("separation = smallest labelled / largest induced X3 difference",
    "(ADEMP at rho = 0: 4.9 for RMTL1, 2.2 for RMTL2)\n\n")
show(separation)

low <- filter(separation, standardised < 2)
if (nrow(low)) {
  cat("WARNING: standardised separation below 2 -",
      paste0(low$target, " at rho = ", low$rho, collapse = "; "), "\n")
}
