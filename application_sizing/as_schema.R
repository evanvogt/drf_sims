##########
# title: application sizing - dataset schemas
##########
# The structure of the two applied datasets the sizing jobs mimic: a platform
# trial analysed as pairwise comparisons within domains (binary and competing
# events outcomes) and a two-arm trial (continuous outcomes). Sourced by
# as_mock.R.
#
# These are for TIMING AND MEMORY ONLY. What drives run time is n, the number
# and type of covariates, how much is missing and the outcome type, so that is
# all the schemas pin down. Arm sizes are rounded to the nearest 50, domains are
# D1-D5 and the event-time profile is pooled over domains. Nothing here is
# meant to reproduce an effect or its heterogeneity.

# ---- platform trial ---------------------------------------------------------
#
# One entry per domain (or domain period). Per domain:
#   arms        n per arm, control first
#   merged      also compare control vs all active arms pooled
#   months      first and last randomisation month (relative); the allocation
#               drifts monthly, so rand_month is the propensity-only covariate
#   n_cont      continuous covariates
#   n_bin       binary covariates
#   cats        one entry per categorical covariate: its number of dummies
#               (missing values are their own category, so never NA)
#   other_doms  one entry per other domain a patient may also be randomised in:
#               its number of arm dummies (all 0 = not randomised there)
#   miss_cont,  how many continuous / binary covariates have any missing
#   miss_bin
#   top_cont,   missing proportion of the single worst continuous / binary
#   top_bin     covariate
#   miss_rho    latent correlation of the missingness indicators within a block
#               (continuous, binary), calibrated by as_mock_check.R so that
#               about CC_TARGET of patients have complete covariates

PLATFORM_DOMAINS <- list(
  D1 = list(arms = c(550, 550, 450), merged = TRUE,  months = c(7, 15),
            n_cont = 28, n_bin = 19, cats = c(2, 3, 4),
            other_doms = c(2, 2, 2, 2, 3),
            miss_cont = 22, miss_bin = 5, top_cont = 0.25, top_bin = 0.25,
            miss_rho = 0.20),
  D2 = list(arms = c(500, 500),      merged = FALSE, months = c(1, 9),
            n_cont = 27, n_bin = 22, cats = c(2, 3, 4),
            other_doms = c(3, 2, 2, 2, 2, 4),
            miss_cont = 21, miss_bin = 7, top_cont = 0.28, top_bin = 0.25,
            miss_rho = 0.22),
  D3 = list(arms = c(250, 250, 200), merged = FALSE, months = c(11, 23),
            n_cont = 30, n_bin = 18, cats = c(2, 3, 4),
            other_doms = c(2, 2, 2),
            miss_cont = 24, miss_bin = 5, top_cont = 0.26, top_bin = 0.22,
            miss_rho = 0.20),
  D4 = list(arms = c(400, 350),      merged = FALSE, months = c(1, 8),
            n_cont = 28, n_bin = 21, cats = c(2, 3, 4),
            other_doms = c(3, 2, 2, 3),
            miss_cont = 23, miss_bin = 5, top_cont = 0.39, top_bin = 0.27,
            miss_rho = 0.27),
  D5 = list(arms = c(300, 600, 450), merged = FALSE, months = c(8, 12),
            n_cont = 28, n_bin = 16, cats = c(2, 3, 4),
            other_doms = c(3, 2, 2, 2),
            miss_cont = 21, miss_bin = 4, top_cont = 0.26, top_bin = 0.28,
            miss_rho = 0.22)
)

# Below the single worst covariate of each type, the missing covariates fall in
# tiers: 6 at 10-20%, 5 at 5-10%, the rest under 5%.
PLATFORM_MISS_TIERS <- list(c(n = 6, lo = 0.10, hi = 0.20),
                            c(n = 5, lo = 0.05, hi = 0.10))
PLATFORM_MISS_LOW <- c(lo = 0.005, hi = 0.045)
CC_TARGET <- 0.22

# Latent correlation of the continuous and binary missingness blocks' factors
MISS_BLOCK_COR <- 0.3

# Competing events outcome, pooled over domains. No censoring: everyone has an
# event by day 90. Cause 1 is discharge, cause 2 death; most deaths are only
# recorded at day 90 (the last follow-up), so no horizon sits at 90.
PLATFORM_HORIZONS <- c(21, 30, 60, 84)
PLATFORM_TIME <- list(
  p_death     = 0.33,                         # D = 2 by day 90
  p_early     = 0.12,                         # share of deaths dated before day 90
  breaks      = c(0, 21, 30, 60, 84, 90),     # day intervals (lower bound open)
  early_death = c(0.80, 0.10, 0.05, 0.05, 0), # deaths dated before day 90, by interval
  discharge   = c(0.79, 0.06, 0.11, 0.025, 0.015) # discharges, by interval
)
# The second binary outcome also counts deaths after discharge
PLATFORM_POST_DC_DEATH <- 0.02

# ---- two-arm trial ----------------------------------------------------------
#
# 1:1 allocation, stratified by a binary site factor (`strata_prob` = 1).
# The baseline value of the outcome (`y_base`, complete) is shared by all
# three outcomes. Outcome means and SDs are rounded to the nearest 5.

TWO_ARM <- list(
  n = 4900, n_cont = 24, n_bin = 7, cats = c(5, 7, 10, 2), strata_prob = 0.25,
  # missing proportions: continuous 3 at 6-8%, 10 at 1-2%, 9 under 1%
  # (2 complete); binary 4 at up to 2% (3 complete)
  miss_cont = c(seq(0.06, 0.08, length.out = 3),
                seq(0.01, 0.02, length.out = 10),
                seq(0.001, 0.009, length.out = 9)),
  miss_bin  = seq(0.005, 0.02, length.out = 4),
  miss_rho  = 0.7,
  outcomes  = data.frame(name = c("y1", "y2", "y3"),
                         mean = c(120, 155, 35),
                         sd   = c(90, 125, 90))
)

# ---- comparisons ------------------------------------------------------------

#' Every pairwise comparison in the platform trial
#'
#' Each active arm against its domain's control, plus the pooled comparison
#' for a domain with `merged = TRUE`.
#'
#' @return data.frame with one row per comparison: id, domain, n_ctrl, n_trt
platform_comparisons <- function(domains = PLATFORM_DOMAINS) {
  rows <- lapply(names(domains), function(d) {
    a <- domains[[d]]$arms
    active <- seq_along(a)[-1]
    out <- data.frame(id = paste0(d, "_", active, "v1"), domain = d,
                      n_ctrl = a[1], n_trt = a[active])
    if (isTRUE(domains[[d]]$merged)) {
      out <- rbind(out, data.frame(
        id = paste0(d, "_", paste(active, collapse = ""), "v1"), domain = d,
        n_ctrl = a[1], n_trt = sum(a[active])))
    }
    out
  })
  out <- do.call(rbind, rows)
  out$n <- out$n_ctrl + out$n_trt
  rownames(out) <- NULL
  out
}

#' Missing proportions of one platform domain's covariates
#'
#' Deterministic (no RNG), so every dataset of a domain shares them.
#'
#' @return list(cont = length n_cont, bin = length n_bin); 0 = never missing
platform_miss_rates <- function(dom) {
  n_rest <- dom$miss_cont + dom$miss_bin - 2
  tiers <- unlist(lapply(PLATFORM_MISS_TIERS, function(t) {
    seq(t[["hi"]], t[["lo"]], length.out = t[["n"]])
  }))
  n_low <- n_rest - length(tiers)
  stopifnot(n_low >= 0)
  rest <- c(tiers, seq(PLATFORM_MISS_LOW[["hi"]], PLATFORM_MISS_LOW[["lo"]],
                       length.out = n_low))
  # the worst of each type first, then the tiers over the continuous
  # covariates, then the binary ones
  cont_rest <- rest[seq_len(dom$miss_cont - 1)]
  bin_rest <- rest[-seq_len(dom$miss_cont - 1)]
  list(
    cont = c(dom$top_cont, cont_rest, rep(0, dom$n_cont - dom$miss_cont)),
    bin  = c(dom$top_bin, bin_rest, rep(0, dom$n_bin - dom$miss_bin))
  )
}
