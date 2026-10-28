##########
# title: csf_sh's wide MCSE under censoring - pre-2026-10-01 tree against now
##########
# In the pre-2026-10-01 results (the DGM before the retune), csf_sh's MSE had
# about twice the pseudo-value CF's MCSE with censoring on and the same MCSE
# without, at about the same mean. The current results do not show it. This
# finds out whether a few tail runs did that, and whether those runs still
# misbehave when csf_sh is refitted on the same saved datasets with today's code:
#
#   1. per run, from the archived tree: MSE of csf_sh and pseudo_cf_whole_oob,
#      the spread of each arm's estimates, and the run's censoring profile
#   2. refit csf_sh on the `--top` censored runs where it does worst against the
#      pseudo-value CF, plus `--top` other runs drawn at random as controls,
#      under get_ipw()'s current rule and its pre-2026-10-10 one (D != 0, which
#      also dropped units censored after the horizon)
#   3. the same summary from the current results (surv_metrics.RDS)
#
# Reading the step 2 output:
#   - each run is refitted `--reps` times with different grf seeds: one csf_sh
#     fit's MSE moves by a few tenths with the seed alone, and the worst runs
#     were picked partly for an unlucky fit, so a single refit would regress
#     to the mean and look like an improvement
#   - archived ~ refit on the controls: the refit reproduces the old fits, so
#     the worst runs can be read off it
#   - worst runs still bad on refit (either rule): it is the data - the old DGM's
#     censored runs - not the code
#   - worst runs fine on refit, and refit_sd_across_seeds large: csf_sh is
#     seed-sensitive on these censored datasets - the tail is the forest's own
#     variance, which the archived single fit happened to land on
#   - worst runs fine on refit with small seed spread: something in the code or
#     the environment (grf 2.4.0 on the cluster) has changed since then
#   - the two rules differ only on the worst runs: the units the old rule dropped
#
# Interactively on the cluster, from competing_risk/, after extracting the
# archive's scenario_<k>/500/censor_<TRUE|FALSE> folders:
#   Rscript csf_sh_tail_check.R --old ~/cr_sept/competing_risk --out ~/csf_sh_tail
# Options: --scenarios 1,3,4  --events 2  --top 10  --reps 3  --threads 1
# Step 2 is 2 x top x reps x 2 rules csf_sh fits per (scenario, event), at a
# few seconds each: about 360 fits at the defaults.
#          --current <surv_metrics.RDS> (default: the live tree's)
# Writes old_per_run.csv, summary_old.csv, refit.csv, summary_refit.csv and
# summary_current.csv to --out.

library(here)
library(dplyr)
library(purrr)
source(here("R", "cate_models.R"))
source(here("competing_risk", "surv_models.R"))
source(here("competing_risk", "surv_config.R"))

args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default) {
  i <- match(paste0("--", name), args)
  if (is.na(i)) default else args[i + 1]
}
old_dir <- opt("old", NA)
if (is.na(old_dir)) stop("--old <extracted pre-2026-10-01 competing_risk/ folder> is required")
old_dir   <- path.expand(old_dir)
out_dir   <- path.expand(opt("out", "csf_sh_tail"))
cur_file  <- path.expand(opt("current", file.path(study$res_path, "surv_metrics.RDS")))
scenarios <- as.integer(strsplit(opt("scenarios", "1,3,4"), ",")[[1]])
events    <- as.integer(strsplit(opt("events", "2"), ",")[[1]])
top       <- as.integer(opt("top", "10"))
reps      <- as.integer(opt("reps", "3"))
threads   <- as.integer(opt("threads", "1"))
horizon   <- 28
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

old_file <- function(scenario, censoring, run) {
  file.path(old_dir, paste0("scenario_", scenario), "500",
            paste0("censor_", censoring), paste0("res_sim_", run, ".RDS"))
}

# share of the spread of x (sum of squared deviations) owed to its k worst runs
top_share <- function(x, k = 5) {
  s <- (x - mean(x))^2
  sum(sort(s, decreasing = TRUE)[seq_len(min(k, length(s)))]) / sum(s)
}

summarise_cells <- function(df) {
  df %>%
    filter(!is.na(csf), !is.na(ps)) %>%
    group_by(across(any_of(c("rho", "scenario", "censoring", "event")))) %>%
    summarise(
      runs       = n(),
      mse_csf    = mean(csf),
      mcse_csf   = sd(csf) / sqrt(n()),
      mse_ps     = mean(ps),
      mcse_ps    = sd(ps) / sqrt(n()),
      mcse_ratio = mcse_csf / mcse_ps,
      diff       = mean(csf - ps),
      mcse_diff  = sd(csf - ps) / sqrt(n()),
      csf_max    = max(csf),
      ps_max     = max(ps),
      csf_top5   = top_share(csf),
      ps_top5    = top_share(ps),
      .groups = "drop"
    )
}

# ---- 1. per run, archived tree ---------------------------------------------
cells <- expand.grid(scenario = scenarios, censoring = c(TRUE, FALSE)) %>%
  mutate(dir = dirname(old_file(scenario, censoring, 1))) %>%
  filter(dir.exists(dir))
if (nrow(cells) == 0) {
  stop("no scenario_<k>/500/censor_<TRUE|FALSE> folders under ", old_dir)
}

read_old_run <- function(f, scenario, censoring) {
  r <- readRDS(f)
  d <- r$data
  run <- as.integer(sub("^res_sim_(\\d+)\\.RDS$", "\\1", basename(f)))
  map_dfr(events, function(ev) {
    est_csf <- r$csf_sh[[paste0("RMST", ev)]]
    est_ps  <- r$pseudo_cf_whole_oob[[paste0("RMTL", ev)]]
    tau_sh  <- r$truth[[paste0("tau_RMST", ev)]]
    tau_ps  <- r$truth[[paste0("tau_RMTL", ev)]]
    tibble(
      scenario, censoring, run, event = ev,
      csf = if (is.null(est_csf)) NA_real_ else mean((est_csf - tau_sh)^2),
      ps  = if (is.null(est_ps)) NA_real_ else mean((est_ps - tau_ps)^2),
      sd_csf = if (is.null(est_csf)) NA_real_ else sd(est_csf),
      sd_ps  = if (is.null(est_ps)) NA_real_ else sd(est_ps),
      sd_tau = sd(tau_ps),
      cens_pre   = sum(d$D == 0 & d$Y < horizon),
      past       = sum(d$Y >= horizon),
      past_trt   = sum(d$Y >= horizon & d$W == 1),
      cens_post  = sum(d$D == 0 & d$Y >= horizon) # dropped by the pre-10-10 rule only
    )
  })
}

message("Step 1: reading ", nrow(cells), " archived cell(s)...")
old_per_run <- pmap_dfr(cells, function(scenario, censoring, dir) {
  files <- list.files(dir, pattern = "^res_sim_\\d+\\.RDS$", full.names = TRUE)
  message("  scenario ", scenario, ", censoring ", censoring, ": ", length(files), " runs")
  map_dfr(files, read_old_run, scenario = scenario, censoring = censoring)
})
summary_old <- summarise_cells(old_per_run)
write.csv(old_per_run, file.path(out_dir, "old_per_run.csv"), row.names = FALSE)
write.csv(summary_old, file.path(out_dir, "summary_old.csv"), row.names = FALSE)
print(summary_old, n = Inf, width = Inf)

# ---- 2. refit csf_sh on the worst and some control censored runs ------------
# csf_sh() looks get_ipw() up in the global environment, so swapping the global
# binding switches the rule; the wrapper also records what the weights did.
get_ipw_current <- get_ipw
ipw_log <- new.env()
use_ipw_rule <- function(rule) {
  assign("get_ipw", function(X, Y, D, W, horizon, censor, num.threads = NULL) {
    w <- get_ipw_current(X, Y, D, W, horizon, censor, num.threads)
    if (rule == "pre_10_10") w$observed <- D != censor
    ipw_log$n_dropped <- sum(!w$observed)
    ipw_log$max_w <- max(w$ipw[w$observed])
    w
  }, envir = globalenv())
}

pick_runs <- function(df) {
  df <- arrange(df, desc(csf - ps))
  worst <- head(df, top)
  rest <- df[-seq_len(nrow(worst)), ]
  control <- rest[sample(nrow(rest), min(top, nrow(rest))), ]
  bind_rows(mutate(worst, group = "worst"), mutate(control, group = "control"))
}

refit_one <- function(p) {
  r <- readRDS(old_file(p$scenario, TRUE, p$run))
  d <- r$data
  X <- as.matrix(d[, !(names(d) %in% c("Y", "D", "W"))])
  tau <- r$truth[[paste0("tau_RMST", p$event)]]
  tidyr::expand_grid(rule = c("current", "pre_10_10"), rep = seq_len(reps)) %>%
    pmap_dfr(function(rule, rep) {
      use_ipw_rule(rule)
      ipw_log$n_dropped <- NA_integer_
      ipw_log$max_w <- NA_real_
      set.seed(100 * p$run + rep) # the same seeds under both rules, so they are paired
      est <- csf_sh(X, d$Y, d$D, d$W, horizon, event = p$event, num.threads = threads)
      tibble(rule, rep, mse_refit = mean((est - tau)^2), sd_refit = sd(est),
             n_dropped = ipw_log$n_dropped, max_w = ipw_log$max_w)
    }) %>%
    bind_cols(select(p, scenario, event, run, group, archived = csf, ps,
                     sd_csf, sd_tau, cens_pre, past, past_trt, cens_post)[rep(1, 2 * reps), ], .)
}

set.seed(1)
picks <- old_per_run %>%
  filter(censoring, !is.na(csf), !is.na(ps)) %>%
  group_split(scenario, event) %>%
  map_dfr(pick_runs)
if (nrow(picks) == 0) {
  message("Step 2: no censored runs in the archive to refit")
  refit <- tibble()
} else {
  message("Step 2: refitting csf_sh on ", nrow(picks), " runs x 2 rules x ", reps, " seeds...")
  refit <- map_dfr(seq_len(nrow(picks)), function(i) {
    p <- picks[i, ]
    message("  scenario ", p$scenario, ", event ", p$event, ", run ", p$run, " (", p$group, ")")
    refit_one(p)
  })
  summary_refit <- refit %>%
    group_by(scenario, event, group, rule) %>%
    summarise(runs = n_distinct(run), archived = mean(archived), refit = mean(mse_refit),
              refit_sd_across_seeds = mean(tapply(mse_refit, run, sd)),
              ps = mean(ps), max_w = median(max_w), n_dropped = mean(n_dropped),
              cens_pre = mean(cens_pre), past = mean(past), .groups = "drop")
  write.csv(refit, file.path(out_dir, "refit.csv"), row.names = FALSE)
  write.csv(summary_refit, file.path(out_dir, "summary_refit.csv"), row.names = FALSE)
  print(summary_refit, n = Inf, width = Inf)
}

# ---- 3. the same summary from the current results ---------------------------
if (!file.exists(cur_file)) {
  message("Step 3: skipped, no ", cur_file)
} else {
  message("Step 3: reading ", cur_file, "...")
  cur <- readRDS(cur_file) %>%
    filter(scenario %in% scenarios, target %in% paste("Event", events)) %>%
    mutate(event = as.integer(sub("Event ", "", target)))
  arm <- function(fw, col) {
    cur %>% filter(framework == fw) %>%
      select(rho, scenario, censoring, run, event, mse) %>%
      rename(all_of(setNames("mse", col)))
  }
  summary_current <- inner_join(arm("csf_sh", "csf"), arm("pseudo_cf_whole_oob", "ps"),
                                by = c("rho", "scenario", "censoring", "run", "event")) %>%
    summarise_cells()
  write.csv(summary_current, file.path(out_dir, "summary_current.csv"), row.names = FALSE)
  print(summary_current, n = Inf, width = Inf)
}

message("Done: ", out_dir)
