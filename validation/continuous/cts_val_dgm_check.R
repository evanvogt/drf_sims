##########
# title: is it the split? - a paired factorial over what changed since August
##########
# A diagnostic, not part of the study pipeline; the companion to
# cts_val_investigation.R. That script compares the archived and current
# results, re-scores the current runs, and tests oracle subgroups without any
# estimator. This one refits the estimators, changing ONE thing at a time
# between the current design and the August one, on the SAME participants.
#
# Why "the same participants" is possible. The rows of a generated trial are
# iid, so the first n1 rows of a trial of n have exactly the distribution of a
# fresh trial of n1 - except for bW, which generate_scenario_data() calibrates
# at the n it is given. The August two-draw design is therefore reproduced
# exactly (in distribution) by splitting one trial and re-calibrating bW inside
# each chunk: Y + W * (bW(n_k) - bW(n)). The other changes are rebuilt the same
# way from the trial's own W, X and noise. So every cell of a run sees the same
# W, X and noise, its estimators are fit from the same seeds, and a difference
# between two cells is that one factor, not seed noise.
#
# Cells (CELLS below), scenario 2 at rho = 0 by default:
#   now        the current design: one trial of 1000 split, bug O baseline,
#              10 covariates, T-learner nuisances
#   chunk_bW   bW calibrated within each chunk - the two-draw design
#   aug_base   the pre-bug-O baseline: Y = 0.3 - 0.05 X1 + 2 X2 + W tau + err,
#              and its bW rule (75% power, sd = s_err + s2). tau's spread is
#              unchanged; X2's prognostic variance goes from 1 to 4
#   aug_covs   August's covariate set: X3 and X5 only where tau uses them (they
#              were drawn only then). In scenario 2, X3 and X5 dropped
#   s_learner  the August S-learner nuisance (nuisance_rf_s below). Only the DR
#              random forest uses it: the causal forest fits its own nuisances,
#              so its rows here must equal `now`'s exactly - a built-in check
#              that the pairing works
#   august     all of the above at once - should reproduce the archived rates
#
# Other scenarios and rho. `scenario=<k>` and `rho=<r>` (a continuous_corr_<r>
# set, so 0 or 0.5) run the same factorial on another CATE, e.g. to see how hard
# a candidate makes subgroup replication before switching the study to it. The
# August cells (aug_base, august) and the comparison with the real runs exist
# only for scenario 2 - the only scenario the study has run, and the one
# baseline AUG below holds - so they are dropped elsewhere; aug_covs is dropped
# where tau uses both X3 and X5, as it would equal `now`. Each scenario and rho
# gets its own folder, cache included. A CATE not in R/dgm_scenarios.R needs a
# row there first: the trial comes from generate_continuous_scenario_data().
#
# The subgroup tests need only chunk 1's tau, so by default only chunk 1 is
# fitted. `oracle_tau` is a third "model": the true tau of chunk 1 through the
# same rpart step, the best any estimator could hand on. Classical and HC3
# p-values are both kept for every test.
#
# `vims` also fits chunk 2 and both importance measures, and scores the
# top-covariate tests with the study's own chunk_validations(). Expect it to
# take more than ten times as long: run it on fewer runs, or cut CELLS down.
#
# Usage, from the repo root or validation/continuous/:
#   Rscript cts_val_dgm_check.R [runs] [workers] [vims] [scenario=<k>] [rho=<r>]
# e.g. Rscript cts_val_dgm_check.R 100 13 scenario=3 rho=0.5
# Defaults: 100 runs, all cores but one, no vims, scenario 2, rho 0. Each run is
# cached under <current metrics folder>/investigation/dgm_check/
# scenario_<k>/rho_<r>/, so an interrupted job resumes, and a finished one
# re-summarises without refitting. Tables print to the console and are written
# there as CSVs.
#
# Run r uses setup_rng_stream(r) and generate_continuous_scenario_data() exactly
# as cts_val_analysis.R does, so in scenario 2 cell `now` is run r of the real
# arm at that rho, participant for participant (its forests are fit from
# different seeds). Every scenario draws the same W, X and noise for a run, so
# the same run number across scenarios is the same participants under a
# different CATE.

library(here)
library(dplyr)
library(tidyr)
library(furrr)
library(grf)
library(rpart)

source(here("R", "utils.R"))
source(here("validation", "continuous", "cts_val_config.R"))
source(here("validation", "continuous", "cts_val_dgms.R"))
source(here("validation", "val_common.R"))

args <- commandArgs(trailingOnly = TRUE)
with_vims <- "vims" %in% args
#' The value of a key=value argument, or `default`
key_arg <- function(key, default) {
  hit <- sub(paste0("^", key, "="), "", grep(paste0("^", key, "="), args, value = TRUE))
  if (length(hit) == 0) default else as.numeric(hit[1])
}
num_args <- suppressWarnings(as.integer(args[args != "vims" & !grepl("=", args)]))
n_runs <- if (length(num_args) >= 1 && !is.na(num_args[1])) num_args[1] else 100L
workers <- if (length(num_args) >= 2 && !is.na(num_args[2])) num_args[2] else
  max(1L, parallel::detectCores() - 1L)

SCENARIO <- key_arg("scenario", 2)
N <- 1000
RHO <- key_arg("rho", 0)
INTERIMS <- c(0.25, 0.5, 0.75)

NOW_PARAMS <- local({
  p <- resolve_set(corr_set("continuous", RHO))
  if (!SCENARIO %in% p$scenario) stop("no scenario ", SCENARIO, " in continuous_corr_", RHO)
  p[p$scenario == SCENARIO, ]
})

# the covariates this scenario's tau uses, read off te_expr - what a variable
# importance measure should put on top
MODIFIERS <- Filter(function(v) te_uses(NOW_PARAMS, v), c("X3", "X4", "X5"))

CELLS <- tribble(
  ~cell,       ~baseline, ~covariates, ~calibration, ~nuisance,
  "now",       "now",     "all",       "trial",      "T",
  "chunk_bW",  "now",     "all",       "chunk",      "T",
  "aug_base",  "aug",     "all",       "trial",      "T",
  "aug_covs",  "now",     "aug",       "trial",      "T",
  "s_learner", "now",     "all",       "trial",      "S",
  "august",    "aug",     "aug",       "chunk",      "S"
)
# see "Other scenarios and rho" in the header
if (SCENARIO != 2) CELLS <- filter(CELLS, baseline != "aug")
if (all(c("X3", "X5") %in% MODIFIERS)) CELLS <- filter(CELLS, covariates != "aug")

# the current metrics and the archived August ones - as cts_val_investigation.R
# finds them
new_candidates <- c(file.path(study$res_path, "cts_val_metrics.RDS"),
                    file.path(dirname(study$res_path), "cts_val_metrics.RDS"))
new_path <- new_candidates[file.exists(new_candidates)][1]
old_tar <- file.path(dirname(here()), "results", "_archive", "pre_2026-09-26",
                     "validation__continuous.tar")
old_member <- "validation/continuous/cts_val_metrics.RDS"

inv_dir <- file.path(if (is.na(new_path)) study$res_path else dirname(new_path),
                     "investigation")
out_dir <- file.path(inv_dir, "dgm_check", paste0("scenario_", SCENARIO),
                     paste0("rho_", RHO))
run_dir <- file.path(out_dir, if (with_vims) "runs_vims" else "runs")
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

#' Print a table and write it to out_dir/<name>.csv
report <- function(df, name, title) {
  cat("\n\n==========", title, "==========\n")
  print(as.data.frame(df), digits = 3, row.names = FALSE)
  write.csv(df, file.path(out_dir, paste0(name, ".csv")), row.names = FALSE)
}

#' Proportion below 0.05 among the non-missing p-values
prop_sig <- function(p) if (all(is.na(p))) NA_real_ else mean(p < 0.05, na.rm = TRUE)

###################
# Building a cell's chunks
###################

# the August baseline of scenario 2 (R/dgm_scenarios.R at d89458b, where it
# was continuous scenario 3). Other scenarios had other baselines then, and
# their August cells are not run
AUG <- list(b0 = 0.3, b1 = -0.05, b2 = 2, s_err = 0.5, s2 = 1)

#' bW under a baseline's own calibration rule, at n
#'
#' August's rule (before bug O): 75% power, sd = s_err + s2, bW = -delta. E[g]
#' is 0 in scenario 2, so the current rule's E[g] adjustment does not enter.
bW_for <- function(baseline, n) {
  if (baseline == "now") return(calibrate_bW(NOW_PARAMS, n, "t"))
  round(-power.t.test(n = n / 2, delta = NULL, sd = AUG$s_err + AUG$s2,
                      power = 0.75)$delta, 2)
}

#' One chunk of the split trial, as a cell would have generated it
#'
#' tau - bW_trial is the bW-free heterogeneity g (-X4 in scenario 2), and
#' Y - p0 - W * tau the trial's own noise, so the outcome is rebuilt from the
#' same participants under the cell's baseline and bW. The current design's
#' chunk is returned untouched rather than rebuilt, so cell `now` is the real
#' study's data to the bit.
#'
#' @param data,truth a chunk from split_trial()
#' @param bW_trial the bW the whole trial was generated with
#' @param cell one row of CELLS
#' @return list(data, tau = the cell's true CATE, bW)
cell_chunk <- function(data, truth, bW_trial, cell) {
  bW <- bW_for(cell$baseline, if (cell$calibration == "chunk") nrow(data) else N)
  tau <- truth$tau - bW_trial + bW

  if (!(cell$baseline == "now" && cell$calibration == "trial")) {
    err <- data$Y - truth$p0 - data$W * truth$tau
    m0 <- if (cell$baseline == "aug") {
      AUG$b0 + AUG$b1 * data$X1 + AUG$b2 * data$X2
    } else {
      truth$p0
    }
    data$Y <- m0 + data$W * tau + err
  }
  if (cell$covariates == "aug") {
    data <- data[, setdiff(names(data), setdiff(c("X3", "X5"), MODIFIERS))]
  }

  list(data = data, tau = tau, bW = bW)
}

###################
# Fitting
###################

#' The August nuisance: one S-learner forest on cbind(W, X)
#'
#' R/cate_models.R's nuisance_rf() as it stood at d89458b, before the move to
#' per-arm forests. Y0.hat / Y1.hat are OOB predictions at the counterfactual
#' W, read through grf's X.orig (grf-labs/grf#307).
nuisance_rf_s <- function(X, Y, W, num.threads = NULL) {
  oob_at <- function(forest, X_counterfactual) {
    forest$X.orig <- X_counterfactual
    forest$predictions <- NULL
    forest$debiased.error <- NULL
    predict(forest)$predictions
  }
  forest <- regression_forest(cbind(W = W, X), Y, num.threads = num.threads)
  Y0.hat <- oob_at(forest, cbind(W = 0, X))
  Y1.hat <- oob_at(forest, cbind(W = 1, X))

  W.hat <- trim_ps(predict(regression_forest(X, W, num.threads = num.threads))$predictions)
  Y.hat.cf <- predict(regression_forest(X, Y, num.threads = num.threads))$predictions

  list(po = dr_pseudo(Y, W, Y1.hat, Y0.hat, W.hat),
       Y.hat = W * Y1.hat + (1 - W) * Y0.hat, Y.hat.cf = Y.hat.cf,
       Y0.hat = Y0.hat, W.hat = W.hat)
}

#' Fit the two forest estimators on one chunk, from a fixed seed
#'
#' The seed is reset before each estimator, so a cell that changes only the
#' nuisance (s_learner) gets the identical causal forest. One grf thread: the
#' parallelism is over runs.
#'
#' @return list(causal_forest, dr_random_forest), each list(tau) - plus te_vims
#'   and shap_vims if vims
fit_chunk <- function(data, nuisance, seed, vims) {
  X <- as.matrix(data[, -c(1, 2)])
  Y <- data$Y
  W <- data$W

  set.seed(seed)
  nuis <- if (nuisance == "S") nuisance_rf_s(X, Y, W, num.threads = 1L) else
    nuisance_rf(X, Y, W, num.threads = 1L)
  set.seed(seed + 1L)
  cf <- run_causal_forest(X, Y, W, nuis, tests = FALSE, num.threads = 1L)
  set.seed(seed + 2L)
  drf <- run_dr_random_forest(X, Y, W, nuis, tests = FALSE, num.threads = 1L)

  fits <- list(causal_forest = list(tau = cf$tau),
               dr_random_forest = list(tau = drf$tau))
  if (vims) {
    set.seed(seed + 3L)
    fits$causal_forest$te_vims <- get_te_vims_causal_forest(X, Y, W, nuis$po, cf$tau,
                                                            num.threads = 1L)
    fits$causal_forest$shap_vims <- get_shap_vims(X, cf$tau)
    set.seed(seed + 4L)
    fits$dr_random_forest$te_vims <- get_te_vims(X, nuis$po, drf$tau, num.threads = 1L)
    fits$dr_random_forest$shap_vims <- get_shap_vims(X, drf$tau)
  }
  fits
}

###################
# Scoring
###################

#' Chunk 1's top/bottom-10% groups, and the tree's prediction of them into chunk 2
#'
#' The subgroup block of chunk_validations() (val_common.R), kept apart here
#' because the diagnostics need the groups themselves, not only the p-values.
#' With vims the p-values are checked against chunk_validations()'s own.
subgroup_split <- function(tau1, X1, X2) {
  rank1 <- rank(tau1, ties.method = "first")
  group <- cut(rank1,
               breaks = quantile(rank1, probs = c(0, 0.1, 0.9, 1)),
               labels = c("bottom10", "middle", "top10"),
               include.lowest = TRUE)
  tree <- rpart(group ~ ., data = data.frame(group = group, X1), method = "class")
  list(group1 = group, pred2 = predict(tree, newdata = X2, type = "class"))
}

#' Subgroup p-values and diagnostics for one model's chunk-1 tau
#'
#' Diagnostics, per side (top / bottom):
#'   cor1        cor(tau_hat, tau) in chunk 1 - how well the estimator ranks
#'   precision1  share of chunk 1's estimated group that is in the true one
#'   n2          chunk-2 participants the tree puts in the group
#'   contrast2   true mean tau inside the group minus outside it in chunk 2 -
#'               the W:v coefficient the test is estimating, i.e. the signal
#'   resid_sd2   residual SD of Y ~ W * v in chunk 2 - the noise
#' Signal over noise is what the rejection rate follows: the estimators move
#' the first, the baseline the second.
score_subgroups <- function(tau_hat1, tau1, c2, X1) {
  data2 <- c2$data
  sp <- subgroup_split(tau_hat1, X1, data2[, -c(1, 2)])
  true_rank1 <- rank(tau1, ties.method = "first")

  pvals <- list()
  diag <- list()
  for (side in c("top", "bottom")) {
    label <- paste0(side, "10")
    v <- as.numeric(sp$pred2 == label)
    true1 <- if (side == "top") true_rank1 > 0.9 * length(tau1) else
      true_rank1 <= 0.1 * length(tau1)
    est1 <- sp$group1 == label
    both <- sum(v) > 0 && sum(v) < length(v)

    pvals[[side]] <- tibble(
      test = paste0("subgroup_", side), se = c("classical", "hc3"),
      p = c(interaction_pval(data2$Y, data2$W, v, robust = FALSE),
            interaction_pval(data2$Y, data2$W, v, robust = TRUE)))
    diag[[side]] <- tibble(
      side = side,
      cor1 = suppressWarnings(cor(tau_hat1, tau1)),
      precision1 = mean(true1[est1]),
      n2 = sum(v),
      contrast2 = if (both) mean(c2$tau[v == 1]) - mean(c2$tau[v == 0]) else NA_real_,
      resid_sd2 = if (both) summary(lm(Y ~ W * v, data.frame(Y = data2$Y, W = data2$W,
                                                             v = v)))$sigma else NA_real_)
  }
  list(pvals = bind_rows(pvals), diag = bind_rows(diag))
}

#' The top-covariate tests from chunk_validations(), long, classical and HC3
#'
#' Also checks chunk_validations()'s subgroup p-values against score_subgroups()'s.
score_vims <- function(fits1, fits2, data1, data2, sub_pvals) {
  out <- list()
  for (robust in c(FALSE, TRUE)) {
    se <- if (robust) "hc3" else "classical"
    val <- chunk_validations(fits1, fits2, data1, data2, robust = robust)
    for (model in names(fits1)) {
      mine <- sub_pvals %>% filter(model == !!model, se == !!se) %>% arrange(test)
      theirs <- unname(val$subgroups[[model]][c("bottom", "top")])
      if (!isTRUE(all.equal(mine$p, theirs))) {
        stop("score_subgroups() disagrees with chunk_validations() for ", model, ", ", se)
      }
      tv <- as.data.frame(val$top_var_tests[[model]])
      # the best-ranked true modifier in each chunk (NA in scenario 1, which has none)
      mod_rank <- val$var_imps[[model]] %>%
        group_by(measure) %>%
        summarise(mod_vi1 = if (any(variables %in% MODIFIERS))
                    max(vi1[variables %in% MODIFIERS]) else NA_real_,
                  mod_vi2 = if (any(variables %in% MODIFIERS))
                    max(vi2[variables %in% MODIFIERS]) else NA_real_,
                  .groups = "drop")
      out[[length(out) + 1]] <- tv %>%
        left_join(mod_rank, by = "measure") %>%
        pivot_longer(c(p_cts, p_cts_adj, p_split), names_to = "p_type", values_to = "p") %>%
        transmute(model = model, test = paste0("topvar_", measure, "_", p_type), se = se, p,
                  x_top, x_top2, p_cov = ncol(data1) - 2, mod_vi1, mod_vi2)
    }
  }
  bind_rows(out)
}

###################
# One run: every interim point and every cell
###################

one_run <- function(r) {
  setup_rng_stream(r)
  gen <- generate_continuous_scenario_data(SCENARIO, N, RHO)

  pvals <- list()
  diag <- list()
  vimrows <- list()
  for (ip in INTERIMS) {
    ch <- split_trial(gen, ip)
    for (k in seq_len(nrow(CELLS))) {
      cell <- CELLS[k, ]
      c1 <- cell_chunk(ch$data1, ch$truth1, gen$bW, cell)
      c2 <- cell_chunk(ch$data2, ch$truth2, gen$bW, cell)
      X1 <- c1$data[, -c(1, 2)]
      keys <- tibble(run = r, interim_prop = ip, cell = cell$cell,
                     bW1 = c1$bW, bW2 = c2$bW)

      # the same seed for this run, interim point and chunk in every cell
      seed1 <- 100000L + 1000L * r + round(100 * ip)
      fits1 <- fit_chunk(c1$data, cell$nuisance, seed1, with_vims)

      taus <- c(lapply(fits1, `[[`, "tau"), list(oracle_tau = c1$tau))
      sub_pvals <- list()
      for (model in names(taus)) {
        s <- score_subgroups(taus[[model]], c1$tau, c2, X1)
        sub_pvals[[model]] <- mutate(s$pvals, model = model)
        diag[[length(diag) + 1]] <- bind_cols(keys, mutate(s$diag, model = model))
      }
      sub_pvals <- bind_rows(sub_pvals)
      pvals[[length(pvals) + 1]] <- bind_cols(keys, sub_pvals)

      if (with_vims) {
        fits2 <- fit_chunk(c2$data, cell$nuisance, seed1 + 200000L, with_vims)
        vimrows[[length(vimrows) + 1]] <- bind_cols(
          keys, score_vims(fits1, fits2, c1$data, c2$data,
                           filter(sub_pvals, model != "oracle_tau")))
      }
    }
  }
  list(cells = CELLS$cell, vims = with_vims, pvals = bind_rows(pvals),
       diag = bind_rows(diag), vimrows = bind_rows(vimrows))
}

run_file <- function(r) file.path(run_dir, paste0("run_", r, ".RDS"))

#' A cached run is reused only if it was made with the same cells
is_cached <- function(r) {
  f <- run_file(r)
  file.exists(f) && identical(readRDS(f)$cells, CELLS$cell)
}

###################
# Fit
###################

todo <- Filter(Negate(is_cached), seq_len(n_runs))
cat(n_runs - length(todo), "of", n_runs, "runs cached;", length(todo), "to fit on",
    workers, "workers", if (with_vims) "(with vims)" else "", "\n")

if (length(todo) > 0) {
  t0 <- Sys.time()
  plan(multisession, workers = workers)
  future_walk(todo, function(r) saveRDS(one_run(r), run_file(r)),
              .options = furrr_options(seed = NULL,
                                       packages = c("dplyr", "tidyr", "grf", "rpart",
                                                    "sandwich", "xgboost",
                                                    "SHAPforxgboost", "furrr")),
              .progress = interactive())
  plan(sequential)
  cat("fitted", length(todo), "runs in",
      format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "\n")
}

###################
# Summarise
###################

runs <- lapply(seq_len(n_runs), function(r) readRDS(run_file(r)))
pvals <- bind_rows(lapply(runs, `[[`, "pvals"))
diag <- bind_rows(lapply(runs, `[[`, "diag"))
vimrows <- bind_rows(lapply(runs, `[[`, "vimrows"))
cell_order <- function(df) mutate(df, cell = factor(cell, levels = CELLS$cell))

report(CELLS %>%
         left_join(pvals %>% distinct(cell, interim_prop, run, bW1, bW2) %>%
                     group_by(cell, interim_prop) %>%
                     summarise(bW1 = first(bW1), bW2 = first(bW2), .groups = "drop"),
                   by = "cell") %>%
         cell_order() %>% arrange(cell, interim_prop),
       "cells", "Cells, and the bW each chunk was generated with")

# ---- rejection rates, one column per cell
rates <- pvals %>%
  group_by(cell, model, test, se, interim_prop) %>%
  summarise(runs = sum(!is.na(p)), reject = prop_sig(p), .groups = "drop")
rates_wide <- rates %>%
  cell_order() %>%
  select(-runs) %>%
  arrange(cell) %>%
  pivot_wider(names_from = cell, values_from = reject) %>%
  arrange(se, test, model, interim_prop)
report(rates_wide, "subgroup_rates",
       paste0("Subgroup tests: proportion p < 0.05 by cell (", n_runs, " runs; MC SE ",
              "about 0.05 at a rate of 0.5)"))

# ---- the decisive table: each cell against `now`, on the same participants
paired <- pvals %>%
  filter(cell != "now") %>%
  inner_join(pvals %>% filter(cell == "now") %>%
               select(run, interim_prop, model, test, se, p_now = p),
             by = c("run", "interim_prop", "model", "test", "se")) %>%
  group_by(cell, model, test, se, interim_prop) %>%
  summarise(pairs = sum(!is.na(p) & !is.na(p_now)),
            reject_now = prop_sig(p_now),
            reject_cell = prop_sig(p),
            diff = reject_cell - reject_now,
            sig_now_only = sum(p_now < 0.05 & p >= 0.05, na.rm = TRUE),
            sig_cell_only = sum(p < 0.05 & p_now >= 0.05, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(mcnemar_p = mapply(function(a, b) {
    if (a + b == 0) NA_real_ else binom.test(b, a + b)$p.value
  }, sig_now_only, sig_cell_only)) %>%
  cell_order() %>%
  arrange(se, test, model, cell, interim_prop)
write.csv(paired, file.path(out_dir, "subgroup_paired_all.csv"), row.names = FALSE)
report(paired %>% filter(se == "hc3"), "subgroup_paired",
       paste("Each cell vs `now` on the same participants (HC3). diff: change in",
             "rejection rate; sig_*_only: runs significant under one only; mcnemar_p:",
             "exact test that diff is 0. s_learner's causal_forest rows must be 0"))

# ---- why: signal (contrast2) and noise (resid_sd2)
diag_sum <- diag %>%
  group_by(cell, model, side, interim_prop) %>%
  summarise(cor1 = mean(cor1, na.rm = TRUE),
            precision1 = mean(precision1, na.rm = TRUE),
            n2 = mean(n2),
            contrast2 = mean(contrast2, na.rm = TRUE),
            resid_sd2 = mean(resid_sd2, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(signal_to_noise = contrast2 / resid_sd2) %>%
  cell_order() %>%
  arrange(side, model, interim_prop, cell)
report(diag_sum, "subgroup_diagnostics",
       paste("Signal and noise behind the subgroup tests: chunk-1 ranking (cor1,",
             "precision1), chunk-2 true contrast (contrast2) and residual SD (resid_sd2)"))

# ---- does the simulation reproduce the real runs at both ends?
#' Real subgroup rejection rates: the arm at RHO (HC3) and August (classical)
#'
#' The arm at RHO comes from the current metrics if they hold it, otherwise
#' from cts_val_investigation.R's subgroups_all.csv in inv_dir, which keeps
#' the rates of whatever metrics it was run against. Scenario 2 only - the
#' study has run no other.
read_real <- function() {
  to_rates <- function(df, source) {
    df %>%
      mutate(interim_prop = as.numeric(as.character(interim_prop))) %>%
      filter(interim_prop %in% INTERIMS) %>%
      pivot_longer(c(top_pval, bottom_pval), names_to = "test", values_to = "p") %>%
      mutate(test = recode(test, top_pval = "subgroup_top",
                           bottom_pval = "subgroup_bottom")) %>%
      group_by(model, test, interim_prop) %>%
      summarise(reject = prop_sig(p), .groups = "drop") %>%
      mutate(source = source)
  }

  real <- list()
  now <- if (!is.na(new_path)) as.data.frame(readRDS(new_path)$subgroups) else NULL
  if (!is.null(now) && "scenario" %in% names(now)) {
    now <- now[as.numeric(as.character(now$scenario)) == SCENARIO, ]
  }
  if (!is.null(now) && "rho" %in% names(now)) {
    rhos <- unique(as.character(now$rho))
    now <- now[as.numeric(as.character(now$rho)) == RHO, ]
  }
  investigation_csv <- file.path(inv_dir, "subgroups_all.csv")
  if (!is.null(now) && nrow(now) > 0) {
    real$now <- to_rates(now, "real_now_hc3")
  } else if (file.exists(investigation_csv)) {
    cat("\nrho =", RHO, "not in the current metrics",
        if (!is.null(now)) paste0("(rho there: ", paste(rhos, collapse = ", "), ")"),
        "- taking it from", investigation_csv, "\n")
    real$now <- read.csv(investigation_csv) %>%
      filter(arm == paste0("now rho=", RHO), interim_prop %in% INTERIMS) %>%
      pivot_longer(c(top_sig, bottom_sig), names_to = "test", values_to = "reject") %>%
      transmute(model, test = recode(test, top_sig = "subgroup_top",
                                     bottom_sig = "subgroup_bottom"),
                interim_prop, reject, source = "real_now_hc3")
  } else {
    cat("\nno real rho =", RHO, "rates found\n")
  }

  if (file.exists(old_tar)) {
    exdir <- tempfile("cts_val_old_")
    dir.create(exdir)
    # tar = "internal": the archive carries bsdtar extended headers that some
    # system tars complain about (and the internal one warns about)
    suppressWarnings(utils::untar(old_tar, files = old_member, exdir = exdir,
                                  tar = "internal"))
    real$aug <- to_rates(as.data.frame(readRDS(file.path(exdir, old_member))$subgroups),
                         "real_august_classical")
  } else {
    cat("\narchived run not found:", old_tar, "\n")
  }
  bind_rows(real)
}

if (SCENARIO == 2) {
  real <- read_real()
  sim_ends <- rates %>%
    filter((cell == "now" & se == "hc3") | (cell == "august" & se == "classical")) %>%
    transmute(source = paste0("sim_", cell, "_", se), model, test, interim_prop, reject)
  report(bind_rows(real, sim_ends) %>%
           pivot_wider(names_from = source, values_from = reject) %>%
           select(model, test, interim_prop, any_of(c("real_now_hc3", "sim_now_hc3",
                                                      "real_august_classical",
                                                      "sim_august_classical"))) %>%
           arrange(test, model, interim_prop),
         "subgroup_vs_real",
         paste0("Calibration: the simulated `now` and `august` cells against the real ",
                "rho = ", RHO, " and archived August rates. If both pairs agree (within ",
                "~0.1 at 100 runs), the cells above account for the whole gap"))
}

# ---- secondary: the top-covariate tests
if (with_vims && nrow(vimrows) > 0) {
  tv <- vimrows %>%
    group_by(cell, model, test, se, interim_prop) %>%
    summarise(reject = prop_sig(p),
              modifier_picked = mean(x_top %in% MODIFIERS),
              modifier_rank1_scaled = mean(mod_vi1 / p_cov),
              chunks_agree = mean(x_top == x_top2),
              .groups = "drop") %>%
    cell_order() %>%
    arrange(se, test, model, interim_prop, cell)
  report(tv %>% filter(se == "hc3"), "topvar_rates",
         paste0("Top-covariate tests (HC3): proportion p < 0.05, how often chunk 1 ",
                "picked a true modifier (", paste(MODIFIERS, collapse = ", "),
                "), the best one's rank / p, and how often the chunks agree on the top"))
  write.csv(tv, file.path(out_dir, "topvar_rates_all.csv"), row.names = FALSE)
}

cat("\nCSVs written to", out_dir, "\n")
