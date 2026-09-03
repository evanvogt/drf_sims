##########
# title: diagnose bias in dr_oracle CATE estimates (continuous study)
##########
# dr_oracle's nuisances (outcome regression, propensity) are the EXACT
# analytical values from the DGP, not estimated - so any bias in its final
# tau should come only from stage 2 (regression_forest smoothing the AIPW
# pseudo-outcome), not from nuisance error. This script decomposes bias into
#   stage1_bias  = mean(po - true_tau)        should be ~0 always
#   final_bias   = mean(tau - true_tau)        production dr_oracle bias
#   stage2_contribution = final_bias - stage1_bias
# and benchmarks stage 2 against a correctly-specified OLS fit of po (which
# should also be ~0-biased, since po is pointwise unbiased for tau) to
# isolate how much of final_bias is regression_forest-specific smoothing.
#
# Run from continuous/:  Rscript cts_oracle_diagnostics.R [n_reps]
# Smoke test first with a small n_reps (e.g. 10-20); the full characterisation
# wants more (e.g. 100) for stable per-cell Monte Carlo error.

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)

source(here("R", "dgm_scenarios.R"))
source(here("R", "cate_models.R"))  # stage2_whole_rf() - the real production stage 2
source(here("R", "metrics.R"))      # cate_metrics() - the real production bias/MSE defs

set <- "continuous"

# ---- config -------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
N_REPS   <- if (length(args) >= 1) as.integer(args[1]) else 20
SCENARIOS <- if (length(args) >= 2) as.integer(strsplit(args[2], ",")[[1]]) else 1:10
NS        <- if (length(args) >= 3) as.integer(strsplit(args[3], ",")[[1]]) else c(100, 250, 500, 1000)

cat(sprintf("dr_oracle diagnostics: %d scenarios x %d sample sizes x %d reps\n",
            length(SCENARIOS), length(NS), N_REPS))

out_dir <- here("continuous", "diagnostics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# correctly-specified stage-2 formula per scenario (RHS only; lm() adds the
# intercept, which is bW). Must match SCENARIO_SETS$continuous$te_expr
# (R/dgm_scenarios.R TE_10) exactly, dropping bW and any b* parameter this
# scenario's te_expr does not actually use (e.g. b5 is a live column in the
# scenario table but TE_10 never references it).
OLS_TERMS <- c(
  "1",                       # 1  no HTE
  "X3",                      # 2
  "X4",                      # 3
  "X3 + X4",                 # 4
  "X3:X4",                   # 5  interaction only, no main effects
  "X3 * X4",                 # 6  main effects + interaction
  "X4:X5",                   # 7  interaction only, no main effects
  "X3 + X4 + X4:X5",         # 8  X3, X4 main effects + X4:X5, no plain X5
  "cos(X4)",                 # 9
  "X3 + I(exp(-abs(X4)))"    # 10
)

# ---- helpers --------------------------------------------------------------

#' Byte-for-byte mirror of run_dr_oracle()'s stage-1 computation
#' (R/cate_models.R:529-546), kept separate only so stage1 (nuisances/po) and
#' stage2 (the forest) can be measured independently. If this ever needs to
#' diverge from run_dr_oracle to keep passing the exact-match check below,
#' that divergence is itself the finding.
#'
#' NOTE: the oracle formula strings reference `X$X1`, `X$X2`, etc. (a literal
#' object named `X`), while the eval envir only flattens X's *columns* into
#' scope - so this lookup resolves only via eval()'s implicit fallback to
#' parent.frame(), where a local object named exactly `X` must be bound. This
#' local's parameter is therefore named `X`, not `Xdf`, to match - renaming it
#' breaks the oracle formula with "object 'X' not found" (confirmed while
#' writing this script), which is itself a live demonstration of how fragile
#' that scoping trick is. See the header comment in R/cate_models.R.
compute_oracle_nuisances <- function(X, W, Y, fmla_info, oracle_link = "identity") {
  link <- if (oracle_link == "logit") plogis else identity
  n_obs <- nrow(X)
  list2env(fmla_info$params, envir = environment())
  fmla <- parse(text = fmla_info$fmla)

  W_temp <- rep(1, n_obs)
  Y1.hat <- link(eval(fmla, envir = list2env(c(list(W = W_temp), X))))
  W_temp <- rep(0, n_obs)
  Y0.hat <- link(eval(fmla, envir = list2env(c(list(W = W_temp), X))))
  Y.hat  <- link(eval(fmla, envir = list2env(c(list(W = W), X))))
  W.hat  <- rep(0.5, n_obs)

  po <- (Y1.hat - Y0.hat) + (Y - Y.hat) * (W - W.hat) / (W.hat * (1 - W.hat))
  list(Y1.hat = Y1.hat, Y0.hat = Y0.hat, Y.hat = Y.hat, W.hat = W.hat, po = po)
}

# ---- one-off sanity check: is the X$X1-via-parent-frame eval trick robust? -
# (see R/cate_models.R header note on run_dr_oracle's formula scoping)
cat("\nChecking oracle formula eval/scoping (scenario 6, n = 500)...\n")
gen_chk <- generate_scenario_data(6, 500, set = set)
fmla_chk <- get_oracle_info(6, gen_chk$bW, set = set)
X_chk <- as.data.frame(as.matrix(gen_chk$dataset[, -c(1, 2)]))
W_chk <- gen_chk$dataset$W
Y_chk <- gen_chk$dataset$Y

nuis_orig <- compute_oracle_nuisances(X_chk, W_chk, Y_chk, fmla_chk)
dev_orig <- max(abs((nuis_orig$Y1.hat - nuis_orig$Y0.hat) - gen_chk$truth$tau))

X_chk_perm <- X_chk[, sample(ncol(X_chk))]  # stress-test: scramble column order
nuis_perm <- compute_oracle_nuisances(X_chk_perm, W_chk, Y_chk, fmla_chk)
dev_perm <- max(abs((nuis_perm$Y1.hat - nuis_perm$Y0.hat) - gen_chk$truth$tau))

cat(sprintf("  max |Y1.hat - Y0.hat - true_tau|: original order = %.3g, permuted order = %.3g\n",
            dev_orig, dev_perm))
if (dev_orig > 1e-8 || dev_perm > 1e-8) {
  warning("Oracle formula does NOT exactly reproduce the true treatment effect - ",
          "this points to a real bug in the oracle nuisance step, not stage-2 smoothing.")
} else {
  cat("  OK - oracle Y1.hat - Y0.hat matches the true treatment effect exactly, ",
      "independent of covariate column order.\n", sep = "")
}

# ---- main simulation loop --------------------------------------------------

cell_results <- list()
cell_i <- 0

for (scenario in SCENARIOS) {
  ols_rhs <- OLS_TERMS[scenario]

  for (n in NS) {
    cat(sprintf("\nscenario %d, n = %d ", scenario, n))

    stage1_bias   <- numeric(N_REPS)
    final_bias    <- numeric(N_REPS)
    ols_bias      <- numeric(N_REPS)
    mse           <- numeric(N_REPS)
    corr          <- numeric(N_REPS)
    te_maxdev     <- numeric(N_REPS)

    for (r in seq_len(N_REPS)) {
      gen <- generate_scenario_data(scenario, n, set = set)
      true_tau <- gen$truth$tau
      fmla_info <- get_oracle_info(scenario, gen$bW, set = set)

      X <- as.matrix(gen$dataset[, -c(1, 2)])
      Xdf <- as.data.frame(X)
      W <- gen$dataset$W
      Y <- gen$dataset$Y

      nuis <- compute_oracle_nuisances(Xdf, W, Y, fmla_info)
      po <- nuis$po

      te_maxdev[r] <- max(abs((nuis$Y1.hat - nuis$Y0.hat) - true_tau))
      stage1_bias[r] <- mean(po - true_tau)

      s <- stage2_whole_rf(X, po)  # real production stage 2
      tau <- s$tau

      m <- cate_metrics(tau, true_tau, scenario)  # real production bias/MSE defs
      final_bias[r] <- m$bias
      mse[r] <- m$mse
      corr[r] <- m$corr

      ols_fit <- lm(as.formula(paste("po ~", ols_rhs)), data = cbind(Xdf, po = po))
      ols_bias[r] <- mean(fitted(ols_fit) - true_tau)

      if (r %% 5 == 0) cat(".")
    }

    cell_i <- cell_i + 1
    cell_results[[cell_i]] <- tibble(
      scenario = scenario, n = n, n_reps = N_REPS,
      stage1_bias_mean = mean(stage1_bias), stage1_bias_mcse = sd(stage1_bias) / sqrt(N_REPS),
      final_bias_mean  = mean(final_bias),  final_bias_mcse  = sd(final_bias) / sqrt(N_REPS),
      ols_bias_mean    = mean(ols_bias),    ols_bias_mcse    = sd(ols_bias) / sqrt(N_REPS),
      stage2_contribution = mean(final_bias) - mean(stage1_bias),
      forest_minus_ols = mean(final_bias) - mean(ols_bias),
      mse_mean = mean(mse), corr_mean = mean(corr),
      te_exact_maxdev = max(te_maxdev)
    )
  }
}
cat("\n")

summary_df <- bind_rows(cell_results)

# ---- interpretation flag ---------------------------------------------------
# Bias not distinguishable from 0 given its MCSE -> fine either way. Bias that
# IS distinguishable from 0 and does not shrink towards the largest n tested
# is flagged for follow-up; shrinking bias (typically scenarios 5-10, small n)
# reads as expected regression-forest smoothing of a nonlinear/interacting
# CATE surface, not a defect.
summary_df <- summary_df %>%
  arrange(scenario, n) %>%
  group_by(scenario) %>%
  mutate(
    z = final_bias_mean / final_bias_mcse,
    abs_bias = abs(final_bias_mean),
    shrank_from_prev = abs_bias < dplyr::lag(abs_bias),
    flag = case_when(
      abs(z) < 2 ~ "not distinguishable from 0",
      n == max(n) & !isTRUE(shrank_from_prev) ~ "PERSISTS at largest n tested - investigate",
      TRUE ~ "shrinking with n - looks like forest smoothing"
    )
  ) %>%
  ungroup() %>%
  select(-z, -abs_bias, -shrank_from_prev)

options(width = 220)
print(as.data.frame(summary_df), row.names = FALSE)

write.csv(summary_df, file.path(out_dir, "cts_oracle_diagnostics_summary.csv"), row.names = FALSE)
saveRDS(summary_df, file.path(out_dir, "cts_oracle_diagnostics_summary.RDS"))

# ---- plot -------------------------------------------------------------------

plot_df <- summary_df %>%
  select(scenario, n, stage1_bias_mean, final_bias_mean, ols_bias_mean) %>%
  pivot_longer(
    cols = c(stage1_bias_mean, final_bias_mean, ols_bias_mean),
    names_to = "component", values_to = "bias"
  ) %>%
  mutate(component = recode(component,
    stage1_bias_mean = "stage 1 (pseudo-outcome)",
    final_bias_mean = "final (dr_oracle, forest stage 2)",
    ols_bias_mean = "OLS benchmark (correct form, stage 2)"
  ))

p <- ggplot(plot_df, aes(x = n, y = bias, color = component)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line() +
  geom_point() +
  facet_wrap(~scenario, scales = "free_y") +
  theme_minimal() +
  labs(title = "dr_oracle CATE bias decomposition, continuous study",
       subtitle = "stage 1 and OLS-benchmark should sit at ~0; the gap to the forest line is stage-2 smoothing",
       x = "n", y = "bias (est - true)", color = NULL)

ggsave(file.path(out_dir, "cts_oracle_bias_decomposition.png"), p, width = 12, height = 7, dpi = 150)

cat(sprintf("\nSaved summary + plot to %s\n", out_dir))
