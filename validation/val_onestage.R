##########
# title: interim-analysis validation - one-stage full-trial comparator, shared
##########
# What does the two-stage approach (nominate a covariate in stage 1, confirm it
# in stage 2) gain over testing every candidate covariate in the whole trial of
# 1000 with a multiplicity correction? This is the one-stage side, run on the
# validation study's own trials, shared by continuous/cts_val_onestage.R and
# binary/bin_val_onestage.R, which call run_onestage(). No CATE model is
# fitted: each trial is OLS only, so the whole thing is one job, not an array.
#
# - Trials: setup_rng_stream(run) then the arm's generator, as
#   <prefix>_val_analysis.R does before split_trial(). The comparator does not
#   depend on interim_prop, so one trial per (scenario, run): scenarios 1-4,
#   runs 1-100, rho and n from the grid.
# - Pairing check: for each trial, the run's results file at interim_prop = 0.5
#   must hold the same trial (rbind of its two chunks). The comparison with the
#   two-stage results is paired by run and rests on this. A mismatch stops the
#   run; a missing file is reported and recorded in `pairing`.
# - Candidates (9): X1-X5, X01-X03, and X04/X05 as one 3-level factor.
# - Tests, all HC3: a joint Wald test of the W: interaction terms of
#   `Y ~ W * g(X_j)`, as an F test on df.residual - for one term that is
#   coef_test()'s HC3 t-test squared (val_common.R), the test p_cts uses. With
#   a binary outcome the model is a linear probability model, linear in the
#   risk, so W: is a risk-difference interaction, as in the two-stage tests.
#   Two forms:
#     spline (primary): g = splines::ns(X_j, df = 4), a fixed-df regression
#       spline, since the shape is estimated in the data it is tested on
#     linear: g = X_j, 1 df, the one-stage p_cts
#   A binary covariate (X1, X3) enters as itself and the factor as its two
#   indicators in both forms.
# - Holm across the 9 p-values per (trial, form), family-wise 0.05. The top
#   covariate is the smallest p, ranked on log p: with strong heterogeneity at
#   n = 1000 several p-values are far below 1e-16, and log p keeps their order
#   even if one would underflow to 0.
#
# run_onestage() writes list(tests, trials, summary):
#   tests    one row per trial x form x covariate
#   trials   one row per trial x form: top covariate and the 0/1 measures
#   summary  per scenario x form x measure: mean over trials and binomial MCSE

library(dplyr)
library(purrr)
library(tibble)
library(sandwich)

ONESTAGE_ALPHA <- 0.05
ONESTAGE_SPLINE_DF <- 4
ONESTAGE_PAIRING_PROP <- 0.5
ONESTAGE_CANDIDATES <- c("X1", "X2", "X3", "X4", "X5", "X01", "X02", "X03",
                         "X04/X05")
ONESTAGE_FORMS <- c("spline", "linear")
ONESTAGE_TRIAL_MEASURES <- c("any_rejected", "top_correct", "any_mod_rejected",
                             "any_nonmod_rejected")

#' The candidate covariates of one trial, the X04/X05 indicators as one factor
candidate_covariates <- function(d) {
  cats <- factor(ifelse(d$X04 == 1, "A", ifelse(d$X05 == 1, "B", "C")),
                 levels = c("C", "A", "B"))
  c(as.list(d[c("X1", "X2", "X3", "X4", "X5", "X01", "X02", "X03")]),
    list("X04/X05" = cats))
}

#' Joint HC3 Wald test of the W: terms of a fitted lm, as an F test
#'
#' @return c(df, stat, p, log_p); stat and p NA if HC3 is not finite (a
#'   leverage-1 row)
wald_hc3 <- function(fit) {
  b <- coef(fit)
  idx <- grep("^W:", names(b))
  no_test <- c(df = length(idx), stat = NA_real_, p = NA_real_, log_p = NA_real_)
  if (length(idx) == 0 || anyNA(b[idx])) return(no_test)
  V <- vcovHC(fit, type = "HC3")[idx, idx, drop = FALSE]
  if (!all(is.finite(V))) return(no_test)
  q <- length(idx)
  stat <- drop(crossprod(b[idx], solve(V, b[idx]))) / q
  c(df = q, stat = stat,
    p = pf(stat, q, fit$df.residual, lower.tail = FALSE),
    log_p = pf(stat, q, fit$df.residual, lower.tail = FALSE, log.p = TRUE))
}

#' One covariate's interaction test in one form
interaction_wald <- function(Y, W, x, form) {
  d <- data.frame(Y = Y, W = W, x = x)
  spline <- form == "spline" && !is.factor(x) && length(unique(x)) > 2
  fit <- if (spline) {
    lm(Y ~ W * splines::ns(x, df = ONESTAGE_SPLINE_DF), data = d)
  } else {
    lm(Y ~ W * x, data = d)
  }
  wald_hc3(fit)
}

#' The run's results file at the pairing interim point
pairing_file <- function(study, scenario, run) {
  combo <- list(rho = unique(study$grid$rho), scenario = scenario,
                n = unique(study$grid$n), interim_prop = ONESTAGE_PAIRING_PROP)
  file.path(combo_dir(study, combo), paste0("res_sim_", run, ".RDS"))
}

#' Is this trial the one the two-stage run split? "match", "missing" or
#' "mismatch"
check_pairing <- function(dataset, file) {
  if (!file.exists(file)) return("missing")
  res <- readRDS(file)
  stacked <- rbind(res$results1$data, res$results2$data)
  rownames(stacked) <- NULL
  if (identical(stacked, dataset)) "match" else "mismatch"
}

#' The true modifiers of each scenario, read off its te_expr as the thesis
#' figures and tables do; none in the null
scenario_modifiers <- function(set, scenarios) {
  sc_table <- resolve_set(set)
  lapply(setNames(scenarios, scenarios), function(s) {
    sc_row <- sc_table[sc_table$scenario == s, ]
    Filter(function(v) te_uses(sc_row, v), paste0("X", 1:5))
  })
}

#' Every test of one trial: one row per form x covariate
onestage_trial <- function(dataset, scenario, run, modifiers, pairing) {
  covs <- candidate_covariates(dataset)
  map_dfr(ONESTAGE_FORMS, function(form) {
    tests <- map_dfr(ONESTAGE_CANDIDATES, function(v) {
      tibble(covariate = v,
             !!!as.list(interaction_wald(dataset$Y, dataset$W, covs[[v]], form)))
    })
    tests %>%
      mutate(p_holm = p.adjust(p, method = "holm",
                               n = length(ONESTAGE_CANDIDATES)),
             rejected = coalesce(p_holm < ONESTAGE_ALPHA, FALSE),
             is_mod = covariate %in% modifiers) %>%
      mutate(scenario = scenario, run = run, form = form, pairing = pairing,
             .before = 1)
  })
}

#' Per-trial measures; NA where undefined (the modifier measures in the null)
onestage_trials <- function(tests) {
  tests %>%
    group_by(scenario, run, form, pairing) %>%
    summarise(
      top_cov = covariate[order(log_p, -stat)[1]],
      top_is_mod = is_mod[order(log_p, -stat)[1]],
      top_rejected = rejected[order(log_p, -stat)[1]],
      any_rejected = any(rejected),
      any_mod_rejected = if (any(is_mod)) any(rejected & is_mod) else NA,
      any_nonmod_rejected = if (any(is_mod)) any(rejected & !is_mod) else NA,
      .groups = "drop"
    ) %>%
    mutate(top_correct = if_else(scenario == 1, NA, top_is_mod & top_rejected))
}

#' Mean over trials and binomial MCSE of each measure, per scenario x form
onestage_summary <- function(tests, trials) {
  binom_summary <- function(d, value) {
    d %>%
      filter(!is.na(.data[[value]])) %>%
      summarise(runs = n(), mean = mean(.data[[value]]),
                mcse = sqrt(mean * (1 - mean) / runs), .groups = "drop")
  }
  summary_trials <- map_dfr(ONESTAGE_TRIAL_MEASURES, function(m) {
    trials %>%
      group_by(scenario, form) %>%
      binom_summary(m) %>%
      mutate(measure = m, covariate = NA_character_)
  })
  # each true modifier rejected, separately (measure 3)
  summary_mods <- tests %>%
    filter(is_mod) %>%
    group_by(scenario, form, covariate) %>%
    binom_summary("rejected") %>%
    mutate(measure = "mod_rejected")

  bind_rows(summary_trials, summary_mods) %>%
    select(scenario, form, measure, covariate, runs, mean, mcse) %>%
    arrange(scenario, form,
            factor(measure, c(ONESTAGE_TRIAL_MEASURES, "mod_rejected")),
            covariate)
}

#' The one-stage comparator over one arm's trials
#'
#' @param study the arm's study_config() (grid, res_path, path columns)
#' @param generate the arm's generator, function(scenario, n, rho)
#' @param set the arm's scenario set, for the true modifiers
#' @param runs the runs to take, 1-100: every scenario has them
#' @param check run the pairing check; FALSE for a local smoke run without the
#'   two-stage results files
#' @param out_file where to save the list(tests, trials, summary)
#' @return that list, invisibly
run_onestage <- function(study, generate, set, runs = 1:100, check = TRUE,
                         out_file) {
  rho <- unique(study$grid$rho)
  n <- unique(study$grid$n)
  stopifnot(length(rho) == 1, length(n) == 1)
  scenarios <- sort(unique(study$grid$scenario))
  modifiers <- scenario_modifiers(set, scenarios)

  trial_keys <- expand.grid(run = runs, scenario = scenarios)
  tests <- map2_dfr(trial_keys$scenario, trial_keys$run, function(s, r) {
    if (r == runs[1]) cat("scenario", s, "\n")
    setup_rng_stream(r)
    dataset <- generate(s, n, rho)$dataset

    pairing <- if (check) check_pairing(dataset, pairing_file(study, s, r)) else
      "not checked"
    if (pairing == "mismatch") {
      stop("scenario ", s, " run ", r, ": the trial drawn here is not the one ",
           "in ", pairing_file(study, s, r), " - the comparison would not be ",
           "paired by run", call. = FALSE)
    }
    onestage_trial(dataset, s, r, modifiers[[as.character(s)]], pairing)
  })

  pairing_counts <- tests %>% distinct(scenario, run, pairing) %>% count(pairing)
  print(pairing_counts)
  if (check && any(pairing_counts$pairing == "missing")) {
    warning("some trials have no two-stage results file at interim_prop = ",
            ONESTAGE_PAIRING_PROP, " to check against - see `pairing`",
            call. = FALSE)
  }

  trials <- onestage_trials(tests)
  summary <- onestage_summary(tests, trials)
  print(as.data.frame(summary), digits = 3)

  onestage <- list(tests = tests, trials = trials, summary = summary)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(onestage, out_file)
  print(paste("one-stage comparator written to", out_file))
  invisible(onestage)
}
