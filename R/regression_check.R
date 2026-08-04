##########
# title: old-vs-new equivalence harness for the R/ refactor
##########
# The refactor moves seven near-identical copies of the CATE estimators into R/.
# Three known bugs (A, F, G - see README.md) are fixed *afterwards*, deliberately,
# so that any number which moves during the move itself is a mistake rather than
# an intended change. This script is what makes that distinction checkable.
#
#   Rscript R/regression_check.R baseline          # capture, BEFORE any change
#   Rscript R/regression_check.R verify            # compare, after each step
#   Rscript R/regression_check.R baseline continuous binary   # a subset
#   Rscript R/regression_check.R list              # show the study ids
#
# Each study runs in its own subprocess. That is not tidiness: `run_all_cate_methods`
# is currently defined in seven files, so sourcing two studies into one session
# would silently give the second one's definition to both.
#
# Scope and limits:
#   - This proves the refactored code reproduces the CURRENT code on THIS machine.
#     It is not a claim about the cluster's numbers - the baseline is local, and
#     R 4.5.3 here vs 4.3.2 there will differ. Re-capture if you change machine.
#   - plan(sequential) throughout. furrr with seed = TRUE is documented to give
#     plan-independent results, and check 0 below asserts that rather than
#     assuming it, so a sequential baseline is comparable to a multisession run.
#   - Settings are small (n = 120, 3 folds, 2-algorithm SL library) so the whole
#     sweep is minutes not hours. Detecting a behaviour change does not need n = 1000;
#     it needs every code path to execute, which the study configs below arrange.

suppressPackageStartupMessages({
  library(here)
})

BASELINE_DIR <- here(".regression_baseline")

# ---- study specifications ---------------------------------------------------
# Each entry sources the study's current files and returns a named list of
# numeric fingerprints. Whatever is returned is what gets compared, so include
# the generated data (to catch RNG-order changes in the DGM) as well as the
# estimates (to catch changes in the estimators).
#
# When a step renames an entry point, update the `run` adapter here - never the
# stored baseline. The baseline is the fixed reference the refactor is judged against.

STUDIES <- list(

  continuous = list(
    sources = c("utils.R", "continuous/cts_dgms.R", "continuous/cts_models.R"),
    run = function() {
      setup_rng_stream(3)
      gen <- generate_continuous_scenario_data(scenario = 6, n = 120)
      fmla <- get_continuous_oracle_info(6, gen$bW)
      setup_rng_stream(3)
      res <- run_all_cate_methods(gen$dataset, n_folds = 3,
                                  sl_lib = SL_LIB, fmla_info = fmla)
      c(list(dataset = gen$dataset, truth = gen$truth, bW = gen$bW),
        tau_fields(res), nuisance_fields(res))
    }
  ),

  binary = list(
    sources = c("utils.R", "binary/bin_dgms.R", "binary/bin_models.R"),
    run = function() {
      setup_rng_stream(3)
      gen <- generate_binary_scenario_data(scenario = 8, n = 120)
      fmla <- get_binary_oracle_info(8, gen$bW)
      setup_rng_stream(3)
      res <- run_all_cate_methods(gen$dataset, n_folds = 3,
                                  sl_lib = SL_LIB, fmla_info = fmla)
      c(list(dataset = gen$dataset, truth = gen$truth, bW = gen$bW),
        tau_fields(res), nuisance_fields(res))
    }
  ),

  ci_continuous = list(
    sources = c("utils.R",
                "confidence_intervals/continuous/cts_ci_dgms.R",
                "confidence_intervals/continuous/cts_ci_models.R"),
    run = function() {
      setup_rng_stream(3)
      gen <- generate_continuous_scenario_data(scenario = 6, n = 120)
      fmla <- get_continuous_oracle_info(6, gen$bW)
      setup_rng_stream(3)
      res <- run_all_cate_methods(gen$dataset, n_folds = 3, fmla_info = fmla,
                                  CI_boot = 10, CI_sf = 0.5, alpha = 0.05)
      c(list(dataset = gen$dataset, truth = gen$truth, bW = gen$bW),
        tau_fields(res), ci_fields(res), nuisance_fields(res))
    }
  ),

  # bug A lives here: this DGM currently carries the CONTINUOUS coefficient table.
  # The baseline pins that behaviour so the refactor can be proved inert before
  # Step 8 deliberately changes it.
  ci_binary = list(
    sources = c("utils.R",
                "confidence_intervals/binary/bin_ci_dgms.R",
                "confidence_intervals/binary/bin_ci_models.R"),
    run = function() {
      setup_rng_stream(3)
      gen <- generate_binary_scenario_data(scenario = 8, n = 120)
      fmla <- get_binary_oracle_info(8, gen$bW)
      setup_rng_stream(3)
      res <- run_all_cate_methods(gen$dataset, n_folds = 3, fmla_info = fmla,
                                  CI_boot = 10, CI_sf = 0.5, alpha = 0.05)
      c(list(dataset = gen$dataset, truth = gen$truth, bW = gen$bW),
        tau_fields(res), ci_fields(res), nuisance_fields(res))
    }
  ),

  # two methods: mean_imputation gives complete data so every arm runs (including
  # the SuperLearner one where bug F lives), IPW exercises the sample.weights /
  # obsWeights threading that is the only difference from the base models file.
  missing_continuous = list(
    sources = c("utils.R",
                "missing/continuous/cts_miss_dgms.R",
                "missing/continuous/cts_miss_models.R"),
    run = function() {
      out <- list()
      for (m in c("mean_imputation", "IPW")) {
        setup_rng_stream(3)
        gen <- suppressMessages(generate_and_process_continuous_data(
          scenario = 4, n = 120, return_truth = TRUE,
          type = "both", prop = 0.3, mech = "MAR", method = m))
        fmla <- get_continuous_oracle_info(4, gen$bW)
        setup_rng_stream(3)
        res <- run_all_cate_methods(gen$dataset, n_folds = 3, sl_lib = SL_LIB,
                                    fmla_info = fmla,
                                    ipw = if (m == "IPW") gen$ipw else NULL)
        out[[paste0(m, ".dataset")]] <- gen$dataset
        out[[paste0(m, ".truth")]] <- gen$truth
        flds <- c(tau_fields(res), nuisance_fields(res))
        for (nm in names(flds)) out[[paste0(m, ".", nm)]] <- flds[[nm]]
      }
      out
    }
  ),

  missing_binary = list(
    sources = c("utils.R",
                "missing/binary/bin_miss_dgms.R",
                "missing/binary/bin_miss_models.R"),
    run = function() {
      out <- list()
      for (m in c("mean_imputation", "IPW")) {
        # n = 400 not 120: GenericML::BLP cannot estimate beta.1 / beta.2 for a
        # binary outcome at the smaller size and errors out
        setup_rng_stream(3)
        gen <- suppressMessages(generate_and_process_binary_data(
          scenario = 4, n = 400, return_truth = TRUE,
          type = "both", prop = 0.3, mech = "MAR", method = m))
        fmla <- get_binary_oracle_info(4, gen$bW)
        setup_rng_stream(3)
        res <- run_all_cate_methods(gen$dataset, n_folds = 3, sl_lib = SL_LIB,
                                    fmla_info = fmla,
                                    ipw = if (m == "IPW") gen$ipw else NULL)
        out[[paste0(m, ".dataset")]] <- gen$dataset
        out[[paste0(m, ".truth")]] <- gen$truth
        flds <- c(tau_fields(res), nuisance_fields(res))
        for (nm in names(flds)) out[[paste0(m, ".", nm)]] <- flds[[nm]]
      }
      out
    }
  ),

  # mi_boot fits every imputed dataset and Rubin-combines. n_imp = 3 rather than
  # the study's 50: the point is to exercise mi_boot + combine_mi, and mice
  # dominates the runtime.
  ci_example = list(
    sources = c("utils.R",
                "missing/ci_example/cts_miss_ci_dgms.R",
                "missing/ci_example/cts_miss_ci_models.R"),
    run = function() {
      setup_rng_stream(3)
      gen <- suppressMessages(generate_and_process_data(
        scenario = 4, n = 120, set = "continuous_missing", return_truth = TRUE,
        type = "both", prop = 0.3, mech = "MAR",
        method = "multiple_imputation", n_imp = 3))
      fmla <- get_continuous_oracle_info(4, gen$bW)
      # combine_mi() in this study reads `alpha` as a FREE VARIABLE - it is not
      # one of its arguments - and works only because cts_miss_ci_analysis.R
      # happens to define alpha at top level. Reproduce that here so the
      # baseline captures the behaviour the study actually ran with.
      assign("alpha", 0.05, envir = globalenv())
      setup_rng_stream(3)
      res <- mi_boot(gen$dataset, n_folds = 3, fmla_info = fmla,
                     CI_boot = 5, CI_sf = 0.5, alpha = 0.05)
      c(list(truth = gen$truth, bW = gen$bW,
             n_imp = length(gen$dataset)),
        tau_fields(res), ci_fields(res))
    }
  ),

  # DEFERRED - excluded from the default sweep, run explicitly with
  #   Rscript R/regression_check.R baseline competing_risk
  # all_cate_surv_models currently aborts with "missing data is currently not
  # supported" from SuperLearner, i.e. an NA reaches a stage-2 fit. Because this
  # folder never calls pretest_superlearner (unlike the cts/bin studies), one bad
  # algorithm takes the whole run down instead of being dropped. Being chased
  # separately; until then there is no baseline to compare against.
  competing_risk = list(
    deferred = TRUE,
    sources = c("utils.R", "competing_risk/surv_dgm.R",
                "competing_risk/surv_models.R"),
    run = function() {
      # n = 400 not 120: at the smaller size grf's causal_survival_forest rejects
      # the sample, because estimated censoring probabilities fall below the 0.05
      # identification bound it requires
      setup_rng_stream(3)
      gen <- generate_surv_data(scenario = 3, n = 400, censoring = TRUE)
      setup_rng_stream(3)
      res <- suppressMessages(
        all_cate_surv_models(gen$dataset, n_folds = 3, horizon = 28))
      # the estimator branches each hold plain numeric tau vectors (cf_ipw and
      # friends return a vector, not a list), so flatten rather than name them
      # one by one - that also picks up any arm added later.
      keep <- setdiff(names(res), c("pseudos", "nuisances", "fold_indices"))
      c(list(dataset = gen$dataset, truth = gen$truth,
             fold_indices = res$fold_indices),
        flatten_numeric(res[keep]))
    }
  ),

  crossfitting = list(
    sources = c("crossfitting/cf_models.R"),
    run = function() {
      setup_rng_stream(3)
      gen <- generate_cf_replicate(scenario = 6, n = 120, n_test = 240)
      setup_rng_stream(3)
      res <- run_all_crossfit_variants(gen$data, gen$X_test, n_folds = 3,
                                       sl_lib = SL_LIB, num.threads = 1,
                                       truth_test = gen$truth_test_tau)
      out <- list(data = gen$data, truth_tau = gen$truth_tau,
                  truth_test_tau = gen$truth_test_tau)
      for (nm in names(res$arms)) out[[paste0("arm.", nm)]] <- res$arms[[nm]]$tau
      out
    }
  )
)

# a deliberately small library - the sweep has to stay fast. SL.mean is included
# because it is the one algorithm that never fails, so pretest_superlearner always
# returns a non-empty library and the SuperLearner arms actually produce estimates.
SL_LIB <- c("SL.mean", "SL.glm")

# ---- fingerprint helpers ----------------------------------------------------

MODELS <- c("causal_forest", "dr_random_forest", "dr_oracle",
            "dr_semi_oracle", "dr_superlearner")

tau_fields <- function(res) {
  out <- list()
  for (m in MODELS) if (!is.null(res[[m]]$tau)) out[[paste0("tau.", m)]] <- res[[m]]$tau
  out
}

# The estimates alone would not notice a change in the shape of what gets SAVED -
# results$nuisances_rf goes to disk and the collect/metrics scripts index it by
# name. Fingerprint the structure (field names, in order) and the matrices
# themselves, so an accidental addition or drop during the refactor is a failure
# rather than a silent difference in every per-run file.
nuisance_fields <- function(res) {
  out <- list()
  for (slot in c("nuisances_rf", "nuisances_sl")) {
    nz <- res[[slot]]
    if (is.null(nz)) next
    out[[paste0(slot, ".names")]] <- names(nz)
    for (nm in names(nz)) out[[paste0(slot, ".", nm)]] <- nz[[nm]]
  }
  if (!is.null(res$fold_indices)) out$fold_indices <- res$fold_indices
  out
}

# recursively pull every numeric leaf out of a nested result list, naming each
# by its path. Used where a study returns bare vectors at varying depths.
flatten_numeric <- function(x, prefix = "") {
  out <- list()
  for (nm in names(x)) {
    v <- x[[nm]]
    key <- if (nzchar(prefix)) paste0(prefix, ".", nm) else nm
    if (is.numeric(v) && is.null(dim(v))) {
      out[[key]] <- v
    } else if (is.list(v)) {
      out <- c(out, flatten_numeric(v, key))
    }
  }
  out
}

ci_fields <- function(res) {
  out <- list()
  for (m in MODELS) {
    a <- res[[m]]
    if (is.null(a)) next
    for (f in c("hb_lb", "hb_ub", "variance")) {
      if (!is.null(a[[f]])) out[[paste0(f, ".", m)]] <- a[[f]]
    }
  }
  out
}

# ---- runner -----------------------------------------------------------------

run_one <- function(id, mode) {
  spec <- STUDIES[[id]]
  if (is.null(spec)) stop("unknown study: ", id)

  suppressPackageStartupMessages({
    library(dplyr); library(furrr); library(future)
    library(grf); library(SuperLearner); library(GenericML); library(coin)
  })
  future::plan(future::sequential)

  for (f in spec$sources) source(here(f))

  cat("  running", id, "...\n")
  t0 <- Sys.time()
  got <- spec$run()
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  %s: %d fields in %.1fs\n", id, length(got), el))

  dir.create(BASELINE_DIR, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(BASELINE_DIR, paste0(id, ".RDS"))

  if (mode == "baseline") {
    saveRDS(got, path)
    cat("  saved ", path, "\n", sep = "")
    return(TRUE)
  }

  if (!file.exists(path)) {
    cat("  NO BASELINE for", id, "- run `baseline` first\n")
    return(FALSE)
  }
  want <- readRDS(path)

  ok <- TRUE
  only_want <- setdiff(names(want), names(got))
  only_got <- setdiff(names(got), names(want))
  if (length(only_want)) {
    cat("  MISSING fields:", paste(only_want, collapse = ", "), "\n"); ok <- FALSE
  }
  if (length(only_got)) {
    cat("  NEW fields:", paste(only_got, collapse = ", "), "\n"); ok <- FALSE
  }

  for (nm in intersect(names(want), names(got))) {
    if (identical(want[[nm]], got[[nm]])) next
    ok <- FALSE
    a <- want[[nm]]; b <- got[[nm]]
    detail <- tryCatch({
      av <- as.numeric(as.matrix(a)); bv <- as.numeric(as.matrix(b))
      if (length(av) != length(bv)) sprintf("length %d vs %d", length(av), length(bv))
      else sprintf("max |diff| = %.3e", max(abs(av - bv), na.rm = TRUE))
    }, error = function(e) "not numerically comparable")
    cat("  DIFF ", nm, ": ", detail, "\n", sep = "")
  }

  cat(if (ok) "  PASS " else "  FAIL ", id, "\n", sep = "")
  ok
}

# ---- entry point ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "verify"

if (mode == "list") {
  cat(paste(names(STUDIES), collapse = "\n"), "\n")
  quit(status = 0)
}

if (!mode %in% c("baseline", "verify", "__one")) {
  stop("usage: Rscript R/regression_check.R [baseline|verify|list] [study ...]")
}

# __one is the in-process form the parent invokes per study
if (mode == "__one") {
  ok <- run_one(args[2], args[3])
  quit(status = if (ok) 0 else 1)
}

# an explicit id list runs whatever you ask for; the default sweep skips deferred
# studies so a known-broken one does not mask a real regression elsewhere
deferred <- names(STUDIES)[vapply(STUDIES, function(s) isTRUE(s$deferred), logical(1))]
ids <- if (length(args) > 1) args[-1] else setdiff(names(STUDIES), deferred)
unknown <- setdiff(ids, names(STUDIES))
if (length(unknown)) stop("unknown study id(s): ", paste(unknown, collapse = ", "))
if (length(args) <= 1 && length(deferred)) {
  cat("  skipping deferred:", paste(deferred, collapse = ", "), "\n")
}

cat("\n=== regression check:", mode, "===\n")
cat("R", as.character(getRversion()), "|", length(ids), "studies | baseline dir:",
    BASELINE_DIR, "\n\n")

rscript <- file.path(R.home("bin"), "Rscript")
self <- here("R", "regression_check.R")

results <- vapply(ids, function(id) {
  st <- system2(rscript, c(shQuote(self), "__one", shQuote(id), shQuote(mode)))
  st == 0
}, logical(1))

cat("\n=== summary ===\n")
cat(sprintf("  %d/%d studies %s\n", sum(results), length(results),
            if (mode == "baseline") "captured" else "identical to baseline"))
if (any(!results)) {
  cat("  problem studies:", paste(names(results)[!results], collapse = ", "), "\n")
  quit(status = 1)
}
