##########
# title: figures for the thesis - DR learners vs causal forest by signal-to-noise
##########
# How the DR learners' CATE accuracy compares with causal forest's as the
# heterogeneity grows relative to the outcome noise, at the primary rho = 0.5
# covariates of the three correlated DGMs. Each point is one run's
#   RMSE(DR) / RMSE(causal forest)
# (below 1, DR did better), on a log scale, against the DGM's
#   SD(tau) / SD(noise).
# The larger outlined point is each group's median.
#
# The pairs, by DR learner family (the facet rows):
#   continuous / binary   DR-RandomForest and DR-SuperLearner vs Causal forest
#                         (sample_size/correlated/, scenarios 1-4 in the main
#                         figure, 5-10 in *_supp - those have 100 runs, not 500)
#   competing risks       DR-RF on PV, DR-SuperLearner on PV and DR-RSF vs
#                         Causal forest on PV (competing_risk/, all 7
#                         scenarios, Event 1 and Event 2): the comparator fits
#                         the same pseudo-values against the same
#                         subdistribution RMTL truth
#
# The x axis comes from the DGM code alone, over the rho = 0.5 population of
# results_processing/dgm_population.R (dgm_covariates()):
#   SD(tau)     continuous / binary: exact, te_moments(). Competing risks: the
#               RMTL CATE of the event over the population's first
#               DGM_N_POP_SURV draws
#   SD(noise)   the irreducible noise given (X, W), averaged over both arms:
#               continuous   s_err
#               binary       sqrt(E[p(1 - p)]), p the risk in either arm. It
#                            moves with n through the calibrated bW, so the
#                            binary x differs slightly by n
#               competing    sqrt(E[Var(L_k | X, W)]), L_k = (h - T)+ 1(D = k)
#               risks        the days lost to event k by the horizon h. No
#                            censoring: both censoring settings share an x
# The four n (or the two censoring settings) would otherwise sit on top of each
# other at one x, so each is drawn shifted by a small fixed offset (X_OFFSET of
# the panel's x range). Each scenario's medians carry its number; in competing
# risks the event is the median's shape, and x is on a log scale (most of its
# scenario x event cells sit below 0.05).
#
# Writes to ../results/thesis_figures/rel_eff/:
#   rel_eff_snr.png        continuous | binary (scenarios 1-4) | competing risks
#   rel_eff_snr_supp.png   continuous | binary, scenarios 5-10
#   rel_eff_snr_{cts,bin}{,_supp}.png, rel_eff_snr_surv.png   one outcome each
# An outcome whose metrics file is missing is skipped with a message, and so is
# any combined figure that needs it.
#
# Local, about a minute (the competing-risk moments are numerical integration):
#   cd results_processing/thesis_figures; Rscript rel_eff_snr.R
# Optional arguments [results root] [figure dir], for a smoke run on mock
# metrics; the results root is the folder holding correlated/ and
# competing_risk/.

library(here)
library(patchwork)
source(here("R", "figures.R"))
source(here("results_processing", "dgm_population.R"))

# paths
args <- commandArgs(trailingOnly = TRUE)
res_root <- if (length(args) >= 1) args[1] else file.path(dirname(here()), "results")
fig_path <- if (length(args) >= 2) {
  args[2]
} else {
  file.path(dirname(here()), "results", "thesis_figures", "rel_eff")
}
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

SCENARIO_SETS_CB <- list(main = 1:4, supp = 5:10)
# the sample-size studies' file prefixes, which also name the figures
PREFIX <- c(continuous = "cts", binary = "bin")

# DR arm -> its family, the facet row. The comparator of each outcome is CF_ARM
DR_FAMILIES <- c(
  dr_random_forest = "RF",
  dr_superlearner = "SL",
  pseudo_dr_whole_oob = "RF",
  sl_dr_whole = "SL",
  rsf_dr_oob = "RSF"
)
CF_ARM <- c(continuous = "causal_forest", binary = "causal_forest",
            surv = "pseudo_cf_whole_oob")
ARM_LABELS <- c(MODEL_LABELS, SURV_ARM_LABELS)

OUTCOME_TITLES <- c(continuous = "Continuous", binary = "Binary",
                    surv = "Competing risks")
X_LAB <- "SD(τ) / SD(noise)"
Y_LAB <- c(
  continuous = "RMSE ratio, DR / causal forest",
  binary = "RMSE ratio, DR / causal forest",
  surv = "RMSE ratio, DR / causal forest on PV"
)
EVENTS <- c(tau_RMTL1 = "Event 1", tau_RMTL2 = "Event 2")

# the n / censoring offset, as a share of the panel's x range
X_OFFSET <- 0.012

# ---- x axis: SD(tau) / SD(noise) ---------------------------------------------

cv <- dgm_covariates()

#' Continuous / binary SD(tau), SD(noise) and their ratio, by scenario and n
cb_snr <- function(cv) {
  bind_rows(lapply(DGM_OUTCOMES, function(outcome) {
    binary <- outcome == "binary"
    set <- corr_set(outcome, DGM_RHO)
    bind_rows(lapply(unlist(SCENARIO_SETS_CB), function(s) {
      p <- dgm_params(outcome, s)
      sd_tau <- sqrt(te_moments(p)$var)
      bind_rows(lapply(DGM_N_LEVELS, function(n) {
        sd_noise <- if (binary) {
          bW <- calibrate_bW(p, n, calibration_for(set))
          tr <- truth_at(p, bW, TRUE, cv$X1, cv$X2, cv$X3, cv$X4, cv$X5)
          sqrt(0.5 * mean(tr$p0 * (1 - tr$p0)) + 0.5 * mean(tr$p1 * (1 - tr$p1)))
        } else {
          p$s_err
        }
        tibble(outcome = outcome, scenario = s, n = n,
               sd_tau = sd_tau, sd_noise = sd_noise)
      }))
    }))
  })) %>%
    mutate(snr = sd_tau / sd_noise)
}

#' The first two moments of L_k = (h - T)+ 1(D = k) given (x, w)
#'
#' E[L_k] = int_0^h (h - t) f_k(t) dt and E[L_k^2] = int_0^h (h - t)^2 f_k(t) dt,
#' f_k = S_all h_k the subdensity of event k (surv_dgm.R's hazards). E[L_k] is
#' the RMTL of event k, which truth_individual() gets by the nested integral of
#' the CIF.
#'
#' @param params one row of survival_scenario_params
#' @param scale1,scale2 the two events' Weibull scales at (x, w)
#' @return c(m1_1, var_1, m1_2, var_2), the mean and variance for each event
l_moments <- function(params, scale1, scale2) {
  h <- params$event_horizon
  shapes <- c(params$shape1, params$shape2)
  scales <- c(scale1, scale2)
  out <- numeric(0)
  for (k in 1:2) {
    f <- function(t) {
      exp(-cum_haz_all(t, params$shape1, params$shape2, scale1, scale2)) *
        hazard(t, shapes[k], scales[k])
    }
    m1 <- integrate(function(t) (h - t) * f(t), 0, h)$value
    m2 <- integrate(function(t) (h - t)^2 * f(t), 0, h)$value
    out[paste0(c("m1_", "var_"), k)] <- c(m1, m2 - m1^2)
  }
  out
}

#' l_moments() in both arms on surv_truth_grid()'s X2 x X1 x X3 grid
#'
#' @return the grid, with per event k: tau_RMTLk (the CATE) and noise_k, the
#'   variance of L_k averaged over the two arms
surv_moment_grid <- function(params, X2_grid = seq(-4.5, 4.5, length.out = 41)) {
  grid <- expand.grid(X2 = X2_grid, X1 = 0:1, X3 = 0:1)
  # the scales as surv_truth() builds them
  log_s1_0 <- log(params$scale1_base) + params$b1_1 * grid$X1 + params$b2_1 * grid$X2
  log_s2_0 <- log(params$scale2_base) + params$b1_2 * grid$X1 + params$b2_2 * grid$X2
  log_s1_1 <- log_s1_0 + params$bW_1 + params$b3_1 * grid$X3
  log_s2_1 <- log_s2_0 + params$bW_2 + params$b3_2 * grid$X3
  arm <- function(ls1, ls2) {
    t(mapply(function(a, b) l_moments(params, exp(a), exp(b)), ls1, ls2))
  }
  m0 <- arm(log_s1_0, log_s2_0)
  m1 <- arm(log_s1_1, log_s2_1)

  # E[L_k] is the RMTL: check it against truth_individual() at a few cells
  for (i in c(1, 62, 123, 164)) {
    ti <- truth_individual(params$shape1, params$shape2,
                           exp(log_s1_1[i]), exp(log_s1_0[i]),
                           exp(log_s2_1[i]), exp(log_s2_0[i]),
                           params$event_horizon)
    stopifnot(
      abs(m1[i, "m1_1"] - ti$RMTL1_1) < 1e-4, abs(m0[i, "m1_1"] - ti$RMTL1_0) < 1e-4,
      abs(m1[i, "m1_2"] - ti$RMTL2_1) < 1e-4, abs(m0[i, "m1_2"] - ti$RMTL2_0) < 1e-4
    )
  }

  mutate(
    grid,
    tau_RMTL1 = m1[, "m1_1"] - m0[, "m1_1"],
    tau_RMTL2 = m1[, "m1_2"] - m0[, "m1_2"],
    noise_1 = 0.5 * (m0[, "var_1"] + m1[, "var_1"]),
    noise_2 = 0.5 * (m0[, "var_2"] + m1[, "var_2"])
  )
}

#' Competing-risk SD(tau), SD(noise) and their ratio, by scenario and event
surv_snr <- function(cv) {
  cv <- cv[seq_len(min(DGM_N_POP_SURV, nrow(cv))), ]
  bind_rows(lapply(survival_scenario_params$scenario, function(s) {
    gt <- surv_moment_grid(
      survival_scenario_params[survival_scenario_params$scenario == s, ]
    )
    bind_rows(lapply(seq_along(EVENTS), function(k) {
      tau <- surv_truth_at_draws(gt, names(EVENTS)[k], cv$X1, cv$X2, cv$X3)
      noise <- surv_truth_at_draws(gt, paste0("noise_", k), cv$X1, cv$X2, cv$X3)
      tibble(outcome = "surv", scenario = s, target = unname(EVENTS[k]),
             sd_tau = sd(tau), sd_noise = sqrt(mean(noise)))
    }))
  })) %>%
    mutate(snr = sd_tau / sd_noise)
}

snr_cb <- cb_snr(cv)
snr_surv <- surv_snr(cv)

message("SD(tau) / SD(noise), continuous and binary:")
print(as.data.frame(snr_cb), digits = 3)
message("SD(tau) / SD(noise), competing risks:")
print(as.data.frame(snr_surv), digits = 3)

# ---- y axis: per-run RMSE ratio ----------------------------------------------

#' readRDS(), or NULL with a message if the file isn't there
read_if_present <- function(file) {
  if (file.exists(file)) return(readRDS(file))
  message("no ", file, " - its figures are skipped")
  NULL
}

#' Each DR arm's per-run RMSE over its comparator's, in the same run
#'
#' @param metrics one row per run and arm; the arm in `arm_col`
#' @param outcome a name of CF_ARM
#' @param keys the columns that identify a run, rho aside
pair_ratios <- function(metrics, outcome, arm_col, keys) {
  m <- filter(metrics, rho == DGM_RHO)
  cf <- m %>%
    filter(.data[[arm_col]] == CF_ARM[[outcome]]) %>%
    select(all_of(keys), rmse_cf = rmse)
  pairs <- m %>%
    filter(.data[[arm_col]] %in% names(DR_FAMILIES)) %>%
    select(all_of(keys), arm = all_of(arm_col), rmse) %>%
    inner_join(cf, by = keys) %>%
    mutate(ratio = rmse / rmse_cf)
  bad <- !is.finite(pairs$ratio) | pairs$ratio <= 0
  if (any(bad)) {
    message(outcome, ": ", sum(bad), " run pairs with an NA or zero RMSE dropped")
  }
  pairs[!bad, ] %>%
    mutate(
      outcome = outcome,
      # facet rows in family order, labelled with the arm's own name
      family = factor(DR_FAMILIES[arm], levels = unique(DR_FAMILIES)),
      model = factor(ARM_LABELS[arm], levels = unique(ARM_LABELS[names(DR_FAMILIES)]))
    )
}

ratios <- list()
for (outcome in DGM_OUTCOMES) {
  m <- read_if_present(file.path(res_root, "correlated", outcome,
                                 paste0(PREFIX[[outcome]], "_corr_metrics.RDS")))
  if (is.null(m)) next
  ratios[[outcome]] <- pair_ratios(m, outcome, "model", c("scenario", "n", "run")) %>%
    left_join(filter(snr_cb, outcome == !!outcome) %>% select(scenario, n, snr),
              by = c("scenario", "n"))
}
m <- read_if_present(file.path(res_root, "competing_risk", "surv_metrics.RDS"))
if (!is.null(m)) {
  m <- filter(m, target %in% EVENTS)
  ratios$surv <- pair_ratios(m, "surv", "framework",
                             c("scenario", "n", "censoring", "target", "run")) %>%
    left_join(select(snr_surv, scenario, target, snr), by = c("scenario", "target"))
}

# ---- figures -----------------------------------------------------------------

#' One outcome's panel: per-run ratios, group medians, one row per DR family
#'
#' @param d the outcome's rows of `ratios`, already filtered to its scenarios
#' @param outcome a name of CF_ARM
snr_plot <- function(d, outcome) {
  surv <- outcome == "surv"
  if (surv) {
    d <- mutate(d,
                group = label_factor(censoring, SURV_CENSORING_LABELS),
                event = factor(target, levels = EVENTS))
    group_lab <- "Censoring"
  } else {
    d <- mutate(d, group = factor(n, levels = DGM_N_LEVELS), event = "all")
    group_lab <- "Sample size"
  }
  # the group's offset from the true x, centred on it. Competing risks' x runs
  # from 0.006 to 0.3 with half its scenario x event cells below 0.05, so it
  # goes on a log scale and the offset is multiplicative. Continuous / binary
  # scenario 1 has x = 0, so they stay linear
  k <- nlevels(d$group)
  shift <- (as.integer(d$group) - (k + 1) / 2) * X_OFFSET
  d$x <- if (surv) {
    d$snr * 10^(shift * diff(range(log10(d$snr))))
  } else {
    d$snr + shift * diff(range(d$snr))
  }

  med <- d %>%
    group_by(model, scenario, event, group, x) %>%
    summarise(ratio = median(ratio), .groups = "drop")
  # one scenario label per scenario (x event), above its highest median
  lab <- med %>%
    group_by(model, scenario, event) %>%
    summarise(x = if (surv) exp(mean(log(x))) else mean(x),
              ratio = max(ratio), .groups = "drop")

  p <- ggplot(d, aes(x, ratio, colour = group)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_point(size = 0.4, alpha = 0.1)
  # the medians: filled circles, or for competing risks a shape per event
  p <- if (surv) {
    p +
      geom_point(data = med, aes(fill = group, shape = event), colour = "black",
                 size = 2) +
      scale_shape_manual(values = c(21, 24)) +
      scale_x_log10() +
      guides(fill = "none",
             shape = guide_legend(override.aes = list(fill = "grey60")))
  } else {
    p +
      geom_point(data = med, aes(fill = group), shape = 21, colour = "black",
                 size = 2, show.legend = FALSE)
  }
  p +
    geom_text(data = lab, aes(label = scenario), colour = "black", size = 2.5,
              vjust = -1) +
    facet_wrap(~ model, ncol = 1) +
    scale_y_log10() +
    drf_scale() +
    scale_fill_paletteer_d("rcartocolor::Safe") +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    labs(title = OUTCOME_TITLES[[outcome]], x = X_LAB, y = Y_LAB[[outcome]],
         colour = group_lab, shape = "Event") +
    drf_theme() +
    theme(legend.position = "bottom")
}

figs <- list()
for (outcome in intersect(DGM_OUTCOMES, names(ratios))) {
  prefix <- PREFIX[[outcome]]
  for (set in names(SCENARIO_SETS_CB)) {
    d <- filter(ratios[[outcome]], scenario %in% SCENARIO_SETS_CB[[set]])
    if (nrow(d) == 0) {
      message(outcome, ": no ", set, " scenario metrics - figure skipped")
      next
    }
    name <- paste0(prefix, if (set == "supp") "_supp")
    figs[[name]] <- snr_plot(d, outcome)
    save_fig(paste0("rel_eff_snr_", name, ".png"), fig_path, width = 15,
             height = 15, plot = figs[[name]])
  }
}
if (!is.null(ratios$surv)) {
  figs$surv <- snr_plot(ratios$surv, "surv")
  save_fig("rel_eff_snr_surv.png", fig_path, width = 15, height = 20,
           plot = figs$surv)
}

# the main figure: continuous and binary over the top two rows (RF, SL), so
# their rows line up with competing risks', which has a third (RSF); the
# legends gather in the cell left below them
if (all(c("cts", "bin", "surv") %in% names(figs))) {
  figs$combined <- figs$cts + figs$bin + figs$surv + guide_area() +
    plot_layout(
      design = c(area(1, 1, 2, 1), area(1, 2, 2, 2), area(1, 3, 3, 3),
                 area(3, 1, 3, 2)),
      guides = "collect"
    )
  save_fig("rel_eff_snr.png", fig_path, width = 33, height = 20,
           plot = figs$combined)
} else {
  message("rel_eff_snr.png skipped: it needs all three outcomes' metrics")
}

if (all(c("cts_supp", "bin_supp") %in% names(figs))) {
  figs$combined_supp <- (figs$cts_supp | figs$bin_supp) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  save_fig("rel_eff_snr_supp.png", fig_path, width = 25, height = 15,
           plot = figs$combined_supp)
} else {
  message("rel_eff_snr_supp.png skipped: it needs both outcomes' scenarios 5-10")
}
