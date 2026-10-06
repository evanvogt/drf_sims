##########
# title: application sizing - mock datasets from the schemas
##########
# mock_platform_domain(d) and mock_two_arm() build datasets with the shape of
# as_schema.R's schemas, platform_comparison() splits a domain into its
# pairwise comparisons, and as_cate_data() / as_surv_data() put a dataset in
# the form cate_methods() / all_cate_surv_models() take. Seed with
# setup_rng_stream() before calling, as the simulation studies do.
#
# Timing only. Covariates come from a Gaussian copula (exchangeable 0.3, a
# third of the continuous ones log-normal); missingness is MCAR, correlated
# within the continuous and the binary block. One continuous covariate (c01)
# modifies the treatment effect so a fit has something to find, but nothing
# here is calibrated to an effect size.

source(here::here("application_sizing", "as_schema.R"))

COV_RHO <- 0.3

# ---- building blocks --------------------------------------------------------

#' Latent N(0, 1) columns with exchangeable correlation rho (one-factor model)
latent_block <- function(n, k, rho, factor = rnorm(n)) {
  sqrt(rho) * factor + sqrt(1 - rho) * matrix(rnorm(n * k), n, k)
}

#' Continuous and binary covariates sharing one latent factor
#'
#' @return list(cont = n x n_cont matrix, bin = n x n_bin matrix,
#'   z = the continuous latent, for the outcome model)
draw_covariates <- function(n, n_cont, n_bin) {
  f <- rnorm(n)
  z <- latent_block(n, n_cont, COV_RHO, f)
  zb <- latent_block(n, n_bin, COV_RHO, f)
  cont <- z
  skewed <- seq_len(n_cont) %% 3 == 0
  cont[, skewed] <- exp(0.75 * z[, skewed])
  prev <- seq(0.05, 0.5, length.out = n_bin)
  bin <- 1L * sweep(zb, 2, qnorm(1 - prev), ">")
  colnames(cont) <- sprintf("c%02d", seq_len(n_cont))
  colnames(bin) <- sprintf("b%02d", seq_len(n_bin))
  list(cont = cont, bin = bin, z = z)
}

#' Dummy columns of categorical covariates (reference level dropped)
#'
#' @param k number of dummies per categorical covariate
#' @param p_ref probability of the reference level
draw_categoricals <- function(n, k, prefix = "cat", p_ref = 0.4) {
  out <- lapply(seq_along(k), function(j) {
    lev <- sample.int(k[j] + 1, n, replace = TRUE,
                      prob = c(p_ref, rep((1 - p_ref) / k[j], k[j])))
    d <- sapply(seq_len(k[j]), function(l) 1L * (lev == l + 1))
    d <- matrix(d, n, k[j])
    colnames(d) <- sprintf("%s%d_%d", prefix, j, seq_len(k[j]))
    d
  })
  do.call(cbind, out)
}

#' Missingness mask, MCAR, correlated within each block
#'
#' @param rates list of per-column missing proportions, one element per block
#' @param rho latent correlation within a block
#' @return logical matrix, TRUE = missing, columns in the order of `rates`
draw_miss_mask <- function(n, rates, rho, block_cor = MISS_BLOCK_COR) {
  shared <- rnorm(n)
  masks <- lapply(rates, function(r) {
    f <- sqrt(block_cor) * shared + sqrt(1 - block_cor) * rnorm(n)
    m <- latent_block(n, length(r), rho, f)
    sweep(m, 2, qnorm(r), "<")
  })
  do.call(cbind, masks)
}

draw_times <- function(n, probs, breaks) {
  bin <- sample.int(length(probs), n, replace = TRUE, prob = probs)
  lo <- breaks[bin] + 1
  hi <- breaks[bin + 1]
  as.integer(lo + floor(runif(n) * (hi - lo + 1)))
}

# ---- platform trial ---------------------------------------------------------

#' One mock domain of the platform trial, every arm
#'
#' Arm W = 0 is the control, 1, 2, ... the active arms in schema order. The
#' allocation starts equal and drifts linearly by month, so that each arm's
#' share over the domain is its schema share - the drift rand_month is there
#' to explain. Missing covariates are imputed at this level (within domain, by
#' arm) before platform_comparison() splits out the pairwise comparisons.
#'
#' @param d a domain name of PLATFORM_DOMAINS, e.g. "D1"
#' @return list(X = covariate data.frame (with NAs), ps_X = data.frame of
#'   propensity-only covariates (rand_month), W = arm 0..K-1, outcomes =
#'   data.frame(bin1, bin2, time, status), info = list(id, n, horizons))
mock_platform_domain <- function(d) {
  dom <- PLATFORM_DOMAINS[[d]]
  if (is.null(dom)) stop("unknown domain: ", d)
  n <- sum(dom$arms)
  k <- length(dom$arms)

  m <- sample(dom$months[1]:dom$months[2], n, replace = TRUE)
  frac <- (m - dom$months[1]) / diff(dom$months)
  share <- dom$arms / n
  probs <- outer(frac, 2 * (share - 1 / k)) + 1 / k   # n x k, rows sum to 1
  W <- apply(probs, 1, function(p) sample.int(k, 1, prob = p)) - 1L

  cov <- draw_covariates(n, dom$n_cont, dom$n_bin)
  cats <- draw_categoricals(n, dom$cats)
  other <- draw_categoricals(n, dom$other_doms, prefix = "od", p_ref = 0.5)

  # outcomes: death risk from a prognostic score; every active arm has the
  # same effect, which c01 modifies
  z <- cov$z
  trt <- as.integer(W > 0)
  prog <- 0.5 * z[, 1] + 0.4 * z[, 2] - 0.3 * z[, 3] + 0.3 * cov$bin[, 1]
  lp <- qlogis(PLATFORM_TIME$p_death) + prog + trt * (-0.1 + 0.3 * z[, 1])
  death <- rbinom(n, 1, plogis(lp)) == 1

  tm <- PLATFORM_TIME
  status <- ifelse(death, 2L, 1L)
  time <- rep(90L, n)
  early <- death & runif(n) < tm$p_early
  time[early] <- draw_times(sum(early), tm$early_death, tm$breaks)
  # discharge: sicker patients (higher score) go home later
  dc <- !death
  u <- pnorm(0.5 * scale(prog)[dc] + sqrt(0.75) * rnorm(sum(dc)))
  bin_dc <- findInterval(u, cumsum(tm$discharge)[-length(tm$discharge)]) + 1
  lo <- tm$breaks[bin_dc] + 1
  hi <- tm$breaks[bin_dc + 1]
  time[dc] <- as.integer(lo + floor(runif(sum(dc)) * (hi - lo + 1)))

  outcomes <- data.frame(
    bin1 = as.integer(death | (dc & runif(n) < PLATFORM_POST_DC_DEATH)),
    bin2 = as.integer(death),
    time = time, status = status
  )

  # missingness last, so the outcome draws do not depend on it
  rates <- platform_miss_rates(dom)
  mask <- draw_miss_mask(n, rates, dom$miss_rho)
  X <- cbind(cov$cont, cov$bin)
  X[mask] <- NA

  list(X = data.frame(X, cats, other), ps_X = data.frame(rand_month = m),
       W = W, outcomes = outcomes,
       info = list(id = d, n = n, horizons = PLATFORM_HORIZONS))
}

#' One pairwise comparison out of a domain (mock or imputed)
#'
#' Keeps the domain's control (W = 0) and the comparison's active arm(s),
#' recoded to W = 1. Every element with one row per patient is subset.
#'
#' @param dom a mock_platform_domain() object, possibly with X replaced by an
#'   imputed version
#' @param id a row id of platform_comparisons(), e.g. "D1_2v1" or "D1_23v1"
platform_comparison <- function(dom, id) {
  cmp <- platform_comparisons()
  cmp <- cmp[cmp$id == id, ]
  if (nrow(cmp) != 1) stop("unknown comparison: ", id)
  if (dom$info$id != cmp$domain) stop(id, " is not in domain ", dom$info$id)
  active <- as.integer(strsplit(sub(".*_(\\d+)v1$", "\\1", id), "")[[1]]) - 1L
  keep <- dom$W == 0 | dom$W %in% active
  list(X = dom$X[keep, , drop = FALSE], ps_X = dom$ps_X[keep, , drop = FALSE],
       W = as.integer(dom$W[keep] > 0), outcomes = dom$outcomes[keep, , drop = FALSE],
       info = list(id = id, n = sum(keep), horizons = dom$info$horizons))
}

#' One mock pairwise comparison, straight from a fresh domain
mock_platform <- function(id) {
  d <- platform_comparisons()$domain[platform_comparisons()$id == id]
  if (length(d) != 1) stop("unknown comparison: ", id)
  platform_comparison(mock_platform_domain(d), id)
}

# ---- two-arm trial ----------------------------------------------------------

#' One mock dataset of the two-arm trial
#'
#' @return list(X = covariate data.frame (with NAs; includes strata and
#'   y_base), W, outcomes = data.frame(y1, y2, y3), info)
mock_two_arm <- function() {
  s <- TWO_ARM
  n <- s$n
  W <- rbinom(n, 1, 0.5)
  strata <- rbinom(n, 1, s$strata_prob)

  cov <- draw_covariates(n, s$n_cont, s$n_bin)
  cats <- draw_categoricals(n, s$cats)
  z <- cov$z
  zb <- 0.5 * z[, 1] + sqrt(0.75) * rnorm(n)   # latent baseline outcome
  prog <- 0.4 * z[, 2] - 0.3 * z[, 3] + 0.3 * cov$bin[, 1] + 0.2 * strata
  outcomes <- as.data.frame(lapply(seq_len(nrow(s$outcomes)), function(k) {
    o <- s$outcomes[k, ]
    o$mean + o$sd * (0.6 * zb + prog + sqrt(0.4) * rnorm(n) +
                       W * (0.05 + 0.1 * z[, 1]))
  }))
  names(outcomes) <- s$outcomes$name
  y_base <- s$outcomes$mean[1] + s$outcomes$sd[1] * zb

  rates <- list(cont = c(s$miss_cont, rep(0, s$n_cont - length(s$miss_cont))),
                bin = c(s$miss_bin, rep(0, s$n_bin - length(s$miss_bin))))
  mask <- draw_miss_mask(n, rates, s$miss_rho)
  X <- cbind(cov$cont, cov$bin)
  X[mask] <- NA

  list(X = data.frame(X, cats, strata = strata, y_base = y_base),
       W = W, outcomes = outcomes, info = list(id = "two_arm", n = n))
}

# ---- shaping for the model functions ----------------------------------------

#' data.frame(Y, W, X...) for cate_methods()
as_cate_data <- function(mock, outcome) {
  data.frame(Y = mock$outcomes[[outcome]], W = mock$W, mock$X)
}

#' data.frame(Y, D, W, X...) for all_cate_surv_models()
as_surv_data <- function(mock) {
  data.frame(Y = mock$outcomes$time, D = mock$outcomes$status, W = mock$W,
             mock$X)
}
