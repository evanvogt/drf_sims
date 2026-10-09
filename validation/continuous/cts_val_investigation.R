##########
# title: why the continuous validation results moved - archived August run vs now
##########
# A diagnostic, not part of the study pipeline. The archived run (results made
# 2026-08-25, under _archive/pre_2026-09-26/) and the current one differ in
# several ways at once:
#   - the bug O baseline (2026-09-26): b0 0.3 -> 0.4, b1 -0.05 -> -0.5,
#     b2 2 -> 1. tau = bW - X4 is unchanged, but X1/X2 are left out of every
#     `Y ~ W * v` interaction test, so their variance is that test's noise
#   - nuisance_rf moved from an S-learner forest to one forest per arm
#     (2026-09-27/28) - the DR random forest's pseudo-outcome, and the
#     causal forest's TE-VIM scoring
#   - every scenario draws X1-X5 (2026-09-27): 10 covariates, not 8
#   - correlated covariates, rho = 0.5 (2026-10-08)
#   - one trial split at the interim instead of two separate draws (2026-10-08)
#   - the DR SuperLearner as a third model, HC3 interaction tests (2026-10-08)
#
# Two parts:
#   metrics  the archived and current cts_val_metrics.RDS side by side, per
#            model and interim point, for each of the four chunk comparisons.
#            A rho = 0 arm shows up as its own column set once it is collected.
#   sim      an oracle check that needs no fitted model: the true top/bottom 10%
#            tau subgroup (what a perfect estimator and tree would hand on), and
#            the X4 proxy X5, interaction-tested at chunk-2 sizes under the
#            August DGM and the current one at rho = 0 and 0.5, with classical
#            and HC3 standard errors. Differences there are the DGM and HC3
#            alone; whatever the metrics move by beyond them is the estimators
#            (T-learner nuisances, the two extra covariates).
#
# Usage, from the repo root or validation/continuous/:
#   Rscript cts_val_investigation.R [all|metrics|sim] [sim_reps]
# Defaults: all, 1000 reps per DGM and chunk size. Tables print to the console
# and are written as CSVs to <current metrics folder>/investigation/.

library(here)
library(dplyr)
library(tidyr)
source(here("validation/continuous/cts_val_config.R"))

args <- commandArgs(trailingOnly = TRUE)
part <- if (length(args) >= 1) args[1] else "all"
sim_reps <- if (length(args) >= 2) as.integer(args[2]) else 1000L
stopifnot(part %in% c("all", "metrics", "sim"))

LANDMARKS <- c(0.25, 0.5, 0.75)

# the current metrics: where cts_val_metrics.R writes them, or one level up
new_candidates <- c(file.path(study$res_path, "cts_val_metrics.RDS"),
                    file.path(dirname(study$res_path), "cts_val_metrics.RDS"))
new_path <- new_candidates[file.exists(new_candidates)][1]

old_tar <- file.path(dirname(here()), "results", "_archive", "pre_2026-09-26",
                     "validation__continuous.tar")
old_member <- "validation/continuous/cts_val_metrics.RDS"

out_dir <- file.path(if (is.na(new_path)) study$res_path else dirname(new_path),
                     "investigation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#' Print a table and write it to out_dir/<name>.csv
report <- function(df, name, title) {
  cat("\n\n==========", title, "==========\n")
  print(as.data.frame(df), digits = 3, row.names = FALSE)
  write.csv(df, file.path(out_dir, paste0(name, ".csv")), row.names = FALSE)
}

#' Proportion below 0.05 among the non-missing p-values
prop_sig <- function(p) if (all(is.na(p))) NA_real_ else mean(p < 0.05, na.rm = TRUE)

###################
# Metrics: archived vs current
###################

#' Read the archived metrics straight out of the tar, into a temp folder
read_old_metrics <- function() {
  exdir <- tempfile("cts_val_old_")
  dir.create(exdir)
  # tar = "internal": the archive carries bsdtar extended headers that some
  # system tars complain about
  utils::untar(old_tar, files = old_member, exdir = exdir, tar = "internal")
  readRDS(file.path(exdir, old_member))
}

#' Give both runs the same columns, and an `arm` naming the run and rho
harmonise <- function(metrics, version) {
  lapply(metrics, function(df) {
    df <- as.data.frame(df)
    if (!"rho" %in% names(df)) df$rho <- NA
    df$interim_prop <- as.numeric(as.character(df$interim_prop))
    df$arm <- if (version == "aug") "aug (indep)" else paste0("now rho=", df$rho)
    df
  })
}

run_metrics_part <- function() {
  if (is.na(new_path)) stop("no current cts_val_metrics.RDS in: ",
                            paste(new_candidates, collapse = ", "))
  if (!file.exists(old_tar)) stop("archived run not found: ", old_tar)
  cat("current metrics:", new_path, "\narchived metrics:", old_tar, "\n")

  old <- harmonise(read_old_metrics(), "aug")
  new <- harmonise(readRDS(new_path), "now")
  both <- function(name) bind_rows(old[[name]], new[[name]])

  subgroups <- both("subgroups")
  variances <- both("variances")
  var_imps  <- both("var_imps")
  top_var   <- both("top_var")
  if (!"measure" %in% names(var_imps)) var_imps$measure <- "tevim"
  if (!"p_cts_adj" %in% names(top_var)) top_var$p_cts_adj <- NA_real_

  # ---- coverage: are both runs complete, and over the same grid?
  coverage <- subgroups %>%
    group_by(arm, model) %>%
    summarise(interim_points = n_distinct(interim_prop),
              runs_min = min(table(interim_prop)),
              runs_max = max(table(interim_prop)),
              .groups = "drop")
  report(coverage, "coverage", "Runs per arm and model")

  # ---- subgroups: how often the chunk-1 responder groups replicate
  sg <- subgroups %>%
    group_by(arm, model, interim_prop) %>%
    summarise(runs = n(),
              top_sig = prop_sig(top_pval),
              bottom_sig = prop_sig(bottom_pval),
              top_na = mean(is.na(top_pval)),
              bottom_na = mean(is.na(bottom_pval)),
              .groups = "drop")
  write.csv(sg, file.path(out_dir, "subgroups_all.csv"), row.names = FALSE)
  report(sg %>% filter(interim_prop %in% LANDMARKS) %>% arrange(model, interim_prop, arm),
         "subgroups", "Subgroups: proportion with p < 0.05 (landmark interim points)")

  # ---- variances: var(tau_hat) against the true Var(tau) = 1
  vr <- variances %>%
    group_by(arm, model, interim_prop) %>%
    summarise(mean_vt1 = mean(vt1), mean_vt2 = mean(vt2),
              median_vt1 = median(vt1), median_vt2 = median(vt2),
              mean_change = mean(var_change), sd_change = sd(var_change),
              .groups = "drop")
  write.csv(vr, file.path(out_dir, "variances_all.csv"), row.names = FALSE)
  report(vr %>% filter(interim_prop %in% LANDMARKS) %>% arrange(model, interim_prop, arm),
         "variances", "Variances: estimated var(tau) per chunk (truth = 1)")

  # ---- variable importance, one row per run first. X4's rank is shown as
  # rank / p because p went from 8 to 10 covariates.
  vi_run <- var_imps %>%
    group_by(arm, model, measure, interim_prop, run) %>%
    summarise(p = n(),
              top1 = variables[which.max(vi1)],
              top2 = variables[which.max(vi2)],
              x4_rank1_scaled = vi1[variables == "X4"][1] / n(),
              spearman = suppressWarnings(cor(vi1, vi2, method = "spearman")),
              .groups = "drop")

  vi <- vi_run %>%
    group_by(arm, model, measure, interim_prop) %>%
    summarise(p = first(p),
              x4_top_stage1 = mean(top1 == "X4"),
              x4_top_both = mean(top1 == "X4" & top2 == "X4"),
              mean_x4_rank_scaled = mean(x4_rank1_scaled, na.rm = TRUE),
              mean_spearman = mean(spearman, na.rm = TRUE),
              .groups = "drop")
  write.csv(vi, file.path(out_dir, "var_imps_all.csv"), row.names = FALSE)
  report(vi %>% filter(interim_prop %in% LANDMARKS) %>%
           arrange(model, measure, interim_prop, arm),
         "var_imps", "Variable importance: X4 on top, X4 rank / p, stage 1 vs 2 rank correlation")

  # when X4 is not chunk 1's winner, which covariate is? Proxies of X4 (X1-X3,
  # X5, X01-X03 at rho = 0.5) versus the independent X04/X05 tells correlation
  # apart from noise
  runners_up <- vi_run %>%
    filter(top1 != "X4") %>%
    count(arm, model, measure, top1, name = "runs") %>%
    group_by(arm, model, measure) %>%
    mutate(share_of_non_x4 = runs / sum(runs)) %>%
    arrange(arm, model, measure, desc(runs)) %>%
    ungroup()
  report(runners_up, "var_imps_non_x4_winners",
         "Variable importance: chunk-1 winner when it is not X4 (all interim points)")

  # ---- top covariate carried into chunk 2, split by whether chunk 1 picked X4.
  # If the current run replicates more often only where x_top != X4, that is
  # the proxy effect p_cts_adj is there to catch.
  tv <- top_var %>%
    mutate(x_top_is_x4 = x_top == "X4") %>%
    group_by(arm, model, measure, interim_prop) %>%
    summarise(runs = n(),
              x4_picked = mean(x_top_is_x4),
              cts_sig = prop_sig(p_cts),
              cts_adj_sig = prop_sig(p_cts_adj),
              split_sig = prop_sig(p_split),
              cts_sig_when_x4 = prop_sig(p_cts[x_top_is_x4]),
              cts_sig_when_not_x4 = prop_sig(p_cts[!x_top_is_x4]),
              cts_adj_sig_when_not_x4 = prop_sig(p_cts_adj[!x_top_is_x4]),
              .groups = "drop")
  write.csv(tv, file.path(out_dir, "top_var_all.csv"), row.names = FALSE)
  report(tv %>% filter(interim_prop %in% LANDMARKS) %>%
           arrange(model, measure, interim_prop, arm),
         "top_var", "Top covariate: proportion with p < 0.05, overall and by whether X4 was picked")
}

###################
# Oracle simulation: the DGM and HC3 without any estimator
###################
# Chunk-2 sizes for interim_prop 0.75 / 0.5 / 0.25. The subgroups are the true
# top/bottom 10% of tau within the chunk - the best any estimator plus rpart
# could hand on - so a change in their rejection rate between DGMs is the
# outcome noise, not the CATE fit.

ORACLE_SIZES <- c(250, 500, 750)

#' The archived run's scenario 3 (now 2), as R/dgm_scenarios.R had it before
#' bug O: b0 = 0.3, b1 = -0.05, b2 = 2, b4 = -1, s_err = 0.5, independent
#' covariates. bW is left at 0: W is in every model, so the W:v coefficient is
#' exactly invariant to it. X5 is drawn here (it was not then) only so the proxy
#' test has a column to read; it is independent of X4, as it would have been.
generate_aug_dgm <- function(m) {
  W <- rbinom(m, 1, 0.5)
  X1 <- rbinom(m, 1, 0.4)
  X2 <- rnorm(m)
  X4 <- rnorm(m)
  X5 <- rnorm(m)
  tau <- -X4
  Y <- 0.3 - 0.05 * X1 + 2 * X2 + W * tau + rnorm(m, 0, 0.5)
  list(dataset = data.frame(Y = Y, W = W, X1 = X1, X2 = X2, X4 = X4, X5 = X5),
       tau = tau)
}

#' The current generator, scenario 2, at a given rho
generate_now_dgm <- function(m, rho) {
  gen <- generate_continuous_scenario_data(2, m, rho)
  list(dataset = gen$dataset, tau = gen$truth$tau)
}

#' One chunk's oracle tests, classical and HC3, via the study's own
#' interaction_pval() / interaction_pval_adj() (validation/val_common.R)
oracle_tests <- function(gen) {
  d <- gen$dataset
  r <- rank(gen$tau, ties.method = "first")
  top <- as.numeric(r > 0.9 * length(r))
  bottom <- as.numeric(r <= 0.1 * length(r))
  X <- d[, -c(1, 2)]

  out <- c()
  for (robust in c(FALSE, TRUE)) {
    se <- if (robust) "hc3" else "classical"
    out[paste0("top10_", se)] <- interaction_pval(d$Y, d$W, top, robust)
    out[paste0("bottom10_", se)] <- interaction_pval(d$Y, d$W, bottom, robust)
    out[paste0("x5_cts_", se)] <- interaction_pval(d$Y, d$W, d$X5, robust)
    out[paste0("x5_cts_adj_", se)] <- interaction_pval_adj(d$Y, d$W, X, "X5", robust)
  }
  out["resid_sd_top10"] <- summary(lm(Y ~ W * v, data.frame(Y = d$Y, W = d$W, v = top)))$sigma
  out
}

run_sim_part <- function() {
  # interaction_pval(), interaction_pval_adj(), and through val_common.R the
  # shared DGM; sourced here so the metrics part does not need the model stack
  source(here("validation", "val_common.R"))
  source(here("validation", "continuous", "cts_val_dgms.R"))

  dgms <- list(
    "aug (indep, old baseline)" = generate_aug_dgm,
    "now rho=0"   = function(m) generate_now_dgm(m, 0),
    "now rho=0.5" = function(m) generate_now_dgm(m, 0.5)
  )

  set.seed(20261009)
  sims <- list()
  for (dgm in names(dgms)) {
    for (m in ORACLE_SIZES) {
      cat("oracle sim:", dgm, "chunk 2 size", m, "\n")
      res <- t(replicate(sim_reps, oracle_tests(dgms[[dgm]](m))))
      sims[[length(sims) + 1]] <- data.frame(dgm = dgm, chunk2_n = m, res,
                                             check.names = FALSE)
    }
  }
  sims <- bind_rows(sims)

  summary_tab <- sims %>%
    pivot_longer(-c(dgm, chunk2_n, resid_sd_top10), names_to = "test", values_to = "p") %>%
    group_by(dgm, chunk2_n, test) %>%
    summarise(reject = prop_sig(p), .groups = "drop") %>%
    pivot_wider(names_from = test, values_from = reject) %>%
    left_join(sims %>% group_by(dgm, chunk2_n) %>%
                summarise(mean_resid_sd = mean(resid_sd_top10), .groups = "drop"),
              by = c("dgm", "chunk2_n")) %>%
    arrange(chunk2_n, dgm)

  report(summary_tab, "oracle_sim",
         paste0("Oracle subgroups and X5 proxy: rejection rate at 0.05 (", sim_reps,
                " reps). X5 is never a modifier, so its rates are false replications."))
}

if (part %in% c("all", "metrics")) run_metrics_part()
if (part %in% c("all", "sim")) run_sim_part()

cat("\nCSVs written to", out_dir, "\n")
