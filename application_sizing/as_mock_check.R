##########
# title: application sizing - check the mock datasets against the schemas
##########
# Run from application_sizing/:  Rscript as_mock_check.R
# Draws one dataset per comparison and prints its size, allocation, complete-
# case share and event profile, and re-derives each domain's miss_rho (the
# latent missingness correlation that gives CC_TARGET complete cases). If a
# derived miss_rho differs from as_schema.R's by more than about 0.05, copy it
# in.

library(here)
source(here("R", "utils.R"))
source(here("application_sizing", "as_mock.R"))

# ---- miss_rho calibration ---------------------------------------------------

cc_share <- function(rates, rho, n = 2e4) {
  setup_rng_stream(1)
  mean(rowSums(draw_miss_mask(n, rates, rho)) == 0)
}

cat("miss_rho for", CC_TARGET, "complete cases\n")
for (d in names(PLATFORM_DOMAINS)) {
  rates <- platform_miss_rates(PLATFORM_DOMAINS[[d]])
  rho <- uniroot(function(r) cc_share(rates, r) - CC_TARGET,
                 c(0.01, 0.95), tol = 0.005)$root
  cat(sprintf("  %s  derived %.2f  schema %.2f  (independent: %.2f CC)\n",
              d, rho, PLATFORM_DOMAINS[[d]]$miss_rho, cc_share(rates, 0)))
}

# ---- one dataset per comparison ---------------------------------------------

cmp <- platform_comparisons()
rows <- lapply(seq_len(nrow(cmp)), function(i) {
  setup_rng_stream(i)
  m <- mock_platform(cmp$id[i])
  o <- m$outcomes
  at <- function(h, s) mean(o$time <= h & o$status == s)
  data.frame(
    id = cmp$id[i], n = m$info$n, p = ncol(m$X), trt = mean(m$W),
    cc = mean(complete.cases(m$X)), bin1 = mean(o$bin1), bin2 = mean(o$bin2),
    dc21 = at(21, 1), dc84 = at(84, 1), death84 = at(84, 2),
    death90 = at(90, 2)
  )
})
cat("\nplatform comparisons\n")
print(format(do.call(rbind, rows), digits = 2), row.names = FALSE)

setup_rng_stream(1)
m <- mock_two_arm()
cat("\ntwo-arm\n")
print(data.frame(n = m$info$n, p = ncol(m$X), trt = mean(m$W),
                 cc = mean(complete.cases(m$X))), digits = 2, row.names = FALSE)
print(sapply(m$outcomes, function(y) round(c(mean = mean(y), sd = sd(y)), 1)))
