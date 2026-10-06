##########
# title: application sizing - time one grid row
##########
# Rscript as_analysis.R <row> [workers] [grf_threads]
#
# Builds the row's mock dataset, runs what the row times (as_config.R) and saves
# the elapsed seconds per fit, plus peak memory, under combo_dir(). Nothing
# estimated is kept - only how long it took and how much memory it needed.
#
# workers / grf_threads default to 4 / 2, as the jobscript asks for 8 cores;
# run the analysis on the laptop with the same pair so the times carry over.

suppressPackageStartupMessages({
  library(dplyr)
  library(furrr)
  library(grf)
  library(GenericML)
  library(SuperLearner)
  library(here)
})

source(here("R", "utils.R"))
source(here("R", "cate_models.R"))
source(here("R", "missingness.R"))
source(here("competing_risk", "surv_models.R"))
source(here("application_sizing", "as_mock.R"))
source(here("application_sizing", "as_config.R"))

args <- as.numeric(commandArgs(trailingOnly = TRUE))
i <- args[1]
workers <- if (length(args) >= 2 && !is.na(args[2])) args[2] else 4
grf_threads <- if (length(args) >= 3 && !is.na(args[3])) args[3] else 2

param <- study$grid[i, ]
print(param)

# AS_SMOKE=1: a wiring check, minutes not hours - 400 patients, 3 folds, a
# two-learner library, a handful of imputations / bootstraps, results to a
# temporary directory so they never pass for real ones
SMOKE <- nzchar(Sys.getenv("AS_SMOKE"))
if (SMOKE) {
  N_FOLDS <- 3
  N_IMP <- 2
  CI <- modifyList(CI, list(boot = 10))
  SF_N_SIM <- 2
  SF_CI_BOOT <- 5
  sl_libraries <- function(n) c("SL.mean", "SL.glm")
  study$res_path <- file.path(tempdir(), "as_smoke")
}

Sys.setenv(OMP_NUM_THREADS = grf_threads)
metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)

# ---- helpers ----------------------------------------------------------------

#' Peak resident memory of this process, MB (Linux; NA elsewhere)
peak_rss_mb <- function() {
  f <- "/proc/self/status"
  if (!file.exists(f)) return(NA_real_)
  line <- grep("^VmHWM:", readLines(f), value = TRUE)
  as.numeric(gsub("[^0-9]", "", line)) / 1024
}

fits <- list()

#' Run expr, record its elapsed seconds as one row of `fits`
time_fit <- function(label, expr) {
  t0 <- proc.time()[["elapsed"]]
  value <- expr
  secs <- proc.time()[["elapsed"]] - t0
  steps <- if (is.list(value) && !is.null(value$timings)) unlist(value$timings) else NULL
  fits[[length(fits) + 1]] <<- data.frame(
    fit = label, seconds = secs,
    steps = if (is.null(steps)) NA_character_ else
      paste(names(steps), round(steps, 1), sep = "=", collapse = ";")
  )
  cat(sprintf("  %-20s %8.1f s\n", label, secs))
  invisible(value)
}

#' Covariates in the shape `handling` names (see as_config.R)
shape_covariates <- function(X, handling) {
  method <- switch(handling, completed = "mean_imputation",
                   indicator = "missing_indicator", NULL)
  if (is.null(method) || !anyNA(X)) return(X)
  handle_missingness(X, method)$data
}

cate_fit <- function(Y, W, X, X_ps = NULL, family = gaussian(), ci = NULL,
                     ipw = NULL) {
  data <- data.frame(Y = Y, W = W, X)
  cate_methods(data, n_folds = N_FOLDS, sl_lib = sl_libraries(nrow(data)),
               family = family, ipw = ipw, ci = ci, profile = "full",
               models = CATE_MODELS_APPLIED, X_ps = X_ps,
               num.threads = grf_threads, verbose_timing = TRUE)
}

# ---- data -------------------------------------------------------------------

platform <- param$dataset != "two_arm"
domain <- if (platform) sub("_.*", "", param$dataset) else NA
# one dataset per domain (and one two-arm dataset), whichever row asks for it
setup_rng_stream(if (platform) match(domain, names(PLATFORM_DOMAINS)) else 99)
mock <- if (platform) mock_platform_domain(domain) else mock_two_arm()
if (SMOKE) {
  keep <- seq_len(400)
  mock$X <- mock$X[keep, ]
  mock$W <- mock$W[keep]
  mock$outcomes <- mock$outcomes[keep, ]
  if (platform) mock$ps_X <- mock$ps_X[keep, , drop = FALSE]
}

# ---- the row ----------------------------------------------------------------

if (param$task == "impute") {
  outs <- mock$outcomes[IMPUTE_OUTCOMES[[if (platform) "platform" else "two_arm"]]]
  # handle_missingness() wants the outcome as `Y` for IPW; the other outcomes
  # ride along as ordinary predictors
  names(outs)[1] <- "Y"
  data <- data.frame(outs[1], W = mock$W, outs[-1], mock$X)
  if (platform) data <- cbind(data, mock$ps_X)
  time_fit(param$handling, handle_missingness(data, param$handling, n_imp = N_IMP))

} else {
  X <- shape_covariates(mock$X, param$handling)
  if (platform) {
    dat <- platform_comparison(c(mock[setdiff(names(mock), "X")], list(X = X)),
                               param$dataset)
  } else {
    dat <- mock
    dat$X <- X
  }
  X_ps <- dat$ps_X       # NULL for the two-arm trial
  o <- dat$outcomes

  if (param$task == "binary") {
    ys <- c(list(bin1 = o$bin1, bin2 = o$bin2),
            setNames(lapply(PLATFORM_HORIZONS, function(h) {
              as.integer(o$time <= h & o$status == 1)
            }), paste0("dc_by_", PLATFORM_HORIZONS)))
    for (nm in names(ys)) {
      time_fit(nm, cate_fit(ys[[nm]], dat$W, dat$X, X_ps, family = binomial()))
    }

  } else if (param$task == "surv") {
    models <- if (param$handling == "none") SURV_MODELS_NA else SURV_MODELS_APPLIED
    sdata <- data.frame(Y = o$time, D = o$status, W = dat$W, dat$X)
    for (h in PLATFORM_HORIZONS) {
      time_fit(paste0("rmtl_dc_", h), all_cate_surv_models(
        sdata, n_folds = N_FOLDS, horizon = h,
        sl_library = sl_libraries(nrow(sdata)), num.threads = grf_threads,
        models = models, estimands = SURV_ESTIMANDS_APPLIED, X_ps = X_ps))
    }

  } else if (param$task == "cts") {
    for (nm in TWO_ARM$outcomes$name) {
      if (param$handling %in% c("complete_cases", "IPW")) {
        h <- handle_missingness(data.frame(Y = o[[nm]], W = dat$W, mock$X),
                                param$handling)
        time_fit(nm, cate_fit(h$data$Y, h$data$W, h$data[, -(1:2)], ipw = h$ipw))
      } else {
        time_fit(nm, cate_fit(o[[nm]], dat$W, dat$X))
      }
    }

  } else if (param$task == "ci") {
    nm <- CI_OUTCOME[[if (platform) "platform" else "two_arm"]]
    fam <- if (platform) binomial() else gaussian()
    time_fit(paste0(nm, "_ci"), cate_fit(o[[nm]], dat$W, dat$X, X_ps,
                                         family = fam, ci = CI))

  } else if (param$task == "sf") {
    nm <- CI_OUTCOME[[if (platform) "platform" else "two_arm"]]
    Xm <- as.matrix(dat$X)
    nuis <- nuisance_rf(Xm, o[[nm]], dat$W, num.threads = grf_threads, X_ps = X_ps)
    tau_hat <- stage2_whole_rf(Xm, nuis$po, num.threads = grf_threads)$tau
    time_fit(paste0(nm, "_sf_", param$sf), find_optimal_sf(
      Xm, o[[nm]], dat$W, nuis, tau_hat, sf_grid = param$sf, n_sim = SF_N_SIM,
      CI_boot = SF_CI_BOOT, alpha = CI$alpha))
  }
}

# ---- save -------------------------------------------------------------------

# each worker's lifetime peak (one element per worker under the default
# chunking), so the job total is roughly main + sum(workers)
worker_mb <- unlist(future_map(seq_len(workers), function(k) {
  f <- "/proc/self/status"
  if (!file.exists(f)) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", grep("^VmHWM:", readLines(f), value = TRUE))) / 1024
}))

results <- list(
  param = param,
  fits = cbind(param[rep(1, length(fits)), c("dataset", "task", "handling", "sf")],
               do.call(rbind, fits), row.names = NULL),
  n = if (param$task == "impute") nrow(mock$X) else nrow(dat$X),
  p = if (param$task == "impute") ncol(mock$X) else ncol(dat$X),
  peak_main_mb = peak_rss_mb(),
  peak_workers_mb = worker_mb,
  workers = workers, grf_threads = grf_threads,
  node = Sys.info()[["nodename"]], r_version = R.version.string
)
print(results$fits)

output_dir <- combo_dir(study, param)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(output_dir, paste0("res_sim_", param$run, ".RDS")))

print("Timing completed!")
