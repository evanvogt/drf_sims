##########
# title: investigation of Alaa & Van Der Schaar 2019 influence function corrected PEHE metric
##########
# Two checks, picked by the first argument:
#
#   Rscript me_infl_investigation.R                         # population (default)
#   Rscript me_infl_investigation.R sims [main|strategies|split]
#
# population: the proxy scores on this study's own DGM (scenarios 1-4,
# continuous_corr_0.5 - me_dgms.R) with ORACLE nuisances: mu0/mu1 are the true
# p0/p1 and pi is the true 0.5, so any gap between a score and true_pehe is
# the score's own bias, not nuisance error. N_BIG is large enough that Monte
# Carlo noise is negligible; bW is calibrated at CALIB_N, a study n, so the
# true ATE is the one the study sees.
#
# sims: the same scores on the completed replicates, from the collected
# <prefix>_all.RDS - no refitting, the scores only need each candidate's saved
# tau and each arm's saved mu0_DR/mu1_DR/pi/phi. Prints, per score, how well
# it ranks the 9 candidates against true PEHE: Spearman, top-1 hit rate and
# regret, as me_results.qmd measures them.
#
# Both save to INVESTIGATION_DIR, results/correlated/model_evaluation/
# investigation/: population.RDS, and <prefix>_infl_check.RDS per tree
# (me_, me_strat_, me_split_) - never into a tree's own _metrics.RDS.
#
# me_metrics.R is deliberately NOT sourced: when <prefix>_all.RDS exists (on
# the cluster) sourcing it recomputes and overwrites <prefix>_metrics.RDS.
# calc_infl_score()/calc_dr_risk() are copied from it below instead.

library(here)
source(here("model_evaluation", "me_config.R")) # study, me_study(), CANDIDATE_MODELS
source(here("R", "metrics.R")) # compute_metrics(), dplyr/tibble/tidyr
source(here("model_evaluation", "me_dgms.R"))
source(here("model_evaluation", "me_nuisance.R")) # calculate_pseudos()

N_BIG <- 2e6
CALIB_N <- 500 # the grid also has 250 and 1000
PI_TRUE <- 0.5

# every result of this script, whichever tree it scores
INVESTIGATION_DIR <- file.path(study$res_path, "investigation")

save_investigation <- function(x, file) {
  dir.create(INVESTIGATION_DIR, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(INVESTIGATION_DIR, file)
  saveRDS(x, path)
  message("saved ", path)
}

# ---- scores -----------------------------------------------------------------

# copied verbatim from me_metrics.R - see the header for why it isn't sourced
calc_infl_score <- function(tau_hat, tau, pi, Y, W) {
  A <- W - pi
  C <- pi * (1 - pi)
  B <- 2 * W * A / C
  diff_tau <- tau - tau_hat

  pehe_infl <- (1 - B) * tau^2 + B * Y * diff_tau - A * diff_tau^2 + tau_hat^2

  mean(pehe_infl)
}

# copied verbatim from me_metrics.R
calc_dr_risk <- function(tau_hat, tau_dr) {
  mean((tau_hat - tau_dr)^2)
}

# one-step corrected PEHE against an arbitrary (possibly biased) surrogate
onestep <- function(th, sur, phi) {
  mean((sur - th)^2 + 2 * (sur - th) * (phi - sur))
}

# influence score metric from the curth and van der schaar 2023 ICML paper, this is how it was written in their appendix but the code of their simulation studies shows them using the same formula as the 2019 paper
calc_infl_score_cvds <- function(tau_hat, tau, pi, Y, W) {
  A <- W - pi
  C <- pi * (1 - pi)
  B <- 2 * W / C
  diff_tau <- tau - tau_hat

  pehe_infl <- (1 - B) * tau^2 + B * Y * diff_tau - A * diff_tau^2 + tau_hat^2

  mean(pehe_infl)
}

# corrected version of the infl score:
calc_infl_score_corrected <- function(tau_hat, tau, pi, Y, W) {
  A <- W - pi
  C <- pi * (1 - pi)
  B <- 2 * W * A / C
  diff_tau <- tau - tau_hat

  pehe_infl <- (1 - B) *
    tau^2 +
    B * Y * diff_tau +
    2 * A * diff_tau^2 -
    2 * tau_hat^2

  mean(pehe_infl)
}

# ---- population check ---------------------------------------------------------

investigate_scenario <- function(scenario) {
  gen <- generate_me_scenario_data(scenario, N_BIG, seed = 1, calib_n = CALIB_N)
  Y <- gen$dataset$Y
  W <- gen$dataset$W
  X2 <- gen$dataset$X2
  tau <- gen$truth$tau

  oracle <- data.frame(
    mu0_DR = gen$truth$p0,
    mu1_DR = gen$truth$p1,
    pi = PI_TRUE
  )
  phi <- calculate_pseudos(oracle, Y, W)$phi

  cands <- list(
    truth = tau,
    ate = rep(mean(tau), N_BIG),
    half = mean(tau) + 0.5 * (tau - mean(tau)),
    mu1 = gen$truth$p1,
    tau_plus_X2 = tau + X2
  )

  message(sprintf(
    "scenario %d: bW = %.3f, sd(tau) = %.3f",
    scenario,
    gen$bW,
    sd(tau)
  ))

  bind_rows(lapply(names(cands), function(nm) {
    th <- cands[[nm]]
    tibble::tibble(
      scenario = scenario,
      candidate = nm,
      true_pehe = mean((th - tau)^2),
      dr = calc_dr_risk(th, phi),
      infl_tau = calc_infl_score(th, tau, PI_TRUE, Y, W),
      infl_phi = calc_infl_score(th, phi, PI_TRUE, Y, W),
      #infl_v2_tau = calc_infl_score_corrected(th, tau, PI_TRUE, Y, W),
      #infl_v2_phi = calc_infl_score_corrected(th, phi, PI_TRUE, Y, W),
      infl_cvds_tau = calc_infl_score_cvds(th, tau, PI_TRUE, Y, W),
      infl_cvds_phi = calc_infl_score_cvds(th, phi, PI_TRUE, Y, W),
      onestep_biased_sur = onestep(th, tau + 0.3 * X2, phi)
    )
  }))
}

run_population <- function() {
  res <- bind_rows(lapply(
    sort(unique(study$grid$scenario)),
    investigate_scenario
  ))
  save_investigation(res, "population.RDS")

  # results
  print(mutate(res, across(where(is.double), ~ round(.x, 3))), n = Inf)

  # table sorted by true PEHE
  print(arrange(res, scenario, true_pehe), n = Inf)

  invisible(res)
}

# ---- check on the completed replicates --------------------------------------

#' Every score for one candidate of one run, long: one row per
#' (pipeline, arm, score). Long rather than me_metrics.R's wide
#' <score>_<arm>_<pipeline> columns because "cv_shared" has an underscore in it.
#'
#' The plug-in surrogate tau_tilde is the arm's own mu1_DR - mu0_DR, the same
#' nuisances phi is built from. With tau_tilde = phi instead, onestep() is
#' exactly calc_dr_risk(), so the plug-in is the only surrogate worth running
#' it with. infl and cvds are run with both surrogates - infl with phi is what
#' me_metrics.R reports.
infl_check_per_model <- function(model_res, true_tau, model, sim_res, keys) {
  tau_hat <- model_res$tau
  Y <- sim_res$data$Y
  W <- sim_res$data$W
  true_pehe <- mean((tau_hat - true_tau)^2)

  bind_rows(lapply(names(sim_res$nuisances), function(pipeline) {
    bind_rows(lapply(names(sim_res$nuisances[[pipeline]]), function(arm) {
      df <- sim_res$nuisances[[pipeline]][[arm]]
      tau_tilde <- df$mu1_DR - df$mu0_DR

      scores <- c(
        infl          = calc_infl_score(tau_hat, df$phi, df$pi, Y, W),
        infl05        = calc_infl_score(tau_hat, df$phi05, 0.5, Y, W),
        infl_plugin   = calc_infl_score(tau_hat, tau_tilde, df$pi, Y, W),
        infl_plugin05 = calc_infl_score(tau_hat, tau_tilde, 0.5, Y, W),
        cvds_plugin   = calc_infl_score_cvds(tau_hat, tau_tilde, df$pi, Y, W),
        cvds_plugin05 = calc_infl_score_cvds(tau_hat, tau_tilde, 0.5, Y, W),
        cvds_phi      = calc_infl_score_cvds(tau_hat, df$phi, df$pi, Y, W),
        cvds_phi05    = calc_infl_score_cvds(tau_hat, df$phi05, 0.5, Y, W),
        onestep       = onestep(tau_hat, tau_tilde, df$phi),
        onestep05     = onestep(tau_hat, tau_tilde, df$phi05),
        dr            = calc_dr_risk(tau_hat, df$phi),
        dr05          = calc_dr_risk(tau_hat, df$phi05)
      )

      tibble(
        pipeline = pipeline,
        arm = arm,
        score = names(scores),
        value = unname(scores),
        true_pehe = true_pehe
      )
    }))
  }))
}

#' Per run and score: does ranking the 9 candidates by the score match
#' ranking them by true PEHE? Lower is better for every score here.
#' Run/score combinations with any NA score are dropped (the known
#' holdout/automl NA phi - see README.md, "Propensity trimming").
selection_performance <- function(scores, path_cols) {
  scores |>
    group_by(across(all_of(c(path_cols, "run"))), pipeline, arm, score) |>
    filter(!anyNA(value), n() == length(CANDIDATE_MODELS)) |>
    summarise(
      spearman = suppressWarnings(cor(value, true_pehe, method = "spearman")),
      top1 = model[which.min(value)] == model[which.min(true_pehe)],
      regret = true_pehe[which.min(value)] - min(true_pehe),
      .groups = "drop"
    )
}

run_sims <- function(which = "main") {
  st <- me_study(which)
  all_path <- file.path(st$res_path, paste0(st$prefix, "_all.RDS"))
  if (!file.exists(all_path)) {
    stop(all_path, " not found - run me_collect.R ", which, " first", call. = FALSE)
  }

  scores <- compute_metrics(
    st,
    readRDS(all_path),
    models = CANDIDATE_MODELS,
    per_model = infl_check_per_model
  )
  perf <- selection_performance(scores, st$path_cols)

  summary <- perf |>
    group_by(across(all_of(st$path_cols)), pipeline, arm, score) |>
    summarise(
      runs = n(),
      spearman = mean(spearman, na.rm = TRUE),
      top1 = mean(top1),
      regret = mean(regret),
      .groups = "drop"
    )

  # onestep and dr differ by a candidate-independent term within a run, so
  # their picks should agree in (almost) every run - a check on the algebra
  same_pick <- perf |>
    filter(score %in% c("onestep", "dr")) |>
    select(-spearman, -top1) |>
    pivot_wider(names_from = score, values_from = regret) |>
    summarise(share_same_regret = mean(abs(onestep - dr) < 1e-10, na.rm = TRUE))
  message("onestep vs dr, share of runs with the same pick: ", round(same_pick$share_same_regret, 4))

  save_investigation(
    list(scores = scores, perf = perf, summary = summary),
    paste0(st$prefix, "_infl_check.RDS")
  )

  # averaged over everything but the score, then the full breakdown
  print(
    perf |>
      group_by(score) |>
      summarise(
        spearman = mean(spearman, na.rm = TRUE),
        top1 = mean(top1),
        regret = mean(regret),
        .groups = "drop"
      ) |>
      arrange(regret),
    n = Inf
  )
  print(
    summary |>
      select(-spearman, -top1, -runs) |>
      pivot_wider(names_from = score, values_from = regret) |>
      mutate(across(where(is.double), ~ round(.x, 3))),
    n = Inf,
    width = Inf
  )

  invisible(list(scores = scores, perf = perf, summary = summary))
}

# ---- run ----------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) == 0) "population" else args[1]
switch(
  mode,
  population = run_population(),
  sims = run_sims(if (length(args) >= 2) args[2] else "main"),
  stop("unknown mode '", mode, "' - use 'population' or 'sims'", call. = FALSE)
)
