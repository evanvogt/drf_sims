##########
# title: back-fill the dr_random_forest HTE tests into completed results
##########
# A one-off repair, not part of the simulation pipeline. `profile = "missing"`
# used to set dr_rf_tests = FALSE (R/cate_models.R), so dr_random_forest alone
# was built inline and carried no BLP or independence test while every other arm
# did. The flag is now TRUE, but ../results/missing/{binary,continuous} holds
# 19,800 finished runs produced under the old value.
#
# Those runs do NOT need redoing. Both tests are deterministic - GenericML::BLP
# is an OLS with a sandwich vcov, coin::independence_test under
# teststat = "quadratic" uses the asymptotic null - and every input survives in
# the saved file (nuisances_rf$W.hat / Y0.hat / po, data, and the arm's own
# tau). So the three fields can be recomputed exactly, with no refit and no RNG,
# which is what this does. The test functions are NOT reimplemented here:
# sourcing R/cate_models.R is what guarantees a patched file and a re-run file
# went through the identical code.
#
# One field is not recoverable. run_dr_random_forest() also returns
# `variance = s$variance` from the stage-2 forest, which the inline branch threw
# away. Patched files therefore have no dr_random_forest$variance while newly
# run ones do. Nothing in these studies reads it (their *_metrics.R call only
# cate_metrics() and hte_test_metrics(); there are no interval metrics here), so
# the asymmetry is recorded rather than papered over - see missing/binary/README.md.
#
# Usage, per study, from its own <prefix>_patch.R wrapper:
#   patch_hte_tests(study, combo_idx = i)            # one parameter combination
#   patch_hte_tests(study, dry_run = TRUE)           # whole tree, writes nothing

library(here)
source(here("R", "pipeline.R"))      # combos, combo_dir, RES_PATTERN, run_numbers
source(here("R", "cate_models.R"))   # run_blp_whole, run_independence_test_whole

# The arm this repairs. Deliberately a constant rather than an argument: every
# other arm already carries its tests, and a patch that could be pointed at an
# arbitrary model is a patch that can silently overwrite good values.
PATCH_MODEL <- "dr_random_forest"

#' Where a study's patch manifest lives
#'
#' Directly under res_path, NOT inside any combo_dir(), so neither
#' check_all.R's count_found() nor get_results() can see it - both scan
#' combo_dir(study, combo) for RES_PATTERN and never look at the top level.
patch_manifest_dir <- function(study) {
  file.path(study$res_path, paste0(study$prefix, "_hte_patch"))
}

#' Decide what one result file needs, without computing anything
#'
#' Split out from the patching itself so dry_run costs a readRDS and nothing
#' more, and so the guards can be read in one place.
#'
#' @return one of the status strings below
patch_status_of <- function(res) {
  # Multiple imputation. That branch of the analysis scripts builds `results`
  # fresh out of combine_mi(), which returns only tau and variance, so no
  # nuisances were ever saved and res$data is a list of 50 imputed datasets
  # rather than one data frame. There is nothing on disk to compute a BLP from.
  # This is a real gap in the study, not a defect in this script - see the
  # multiple-imputation note in missing/README.md.
  if (is.null(res$nuisances_rf) || !is.data.frame(res$data)) {
    return("skipped_mi")
  }
  # The arm should always be present - unlike dr_oracle / dr_semi_oracle /
  # dr_superlearner it is not gated on a complete covariate matrix - but a
  # missing arm must be reported, never invented.
  if (is.null(res[[PATCH_MODEL]]) || is.null(res[[PATCH_MODEL]]$tau)) {
    return("skipped_no_arm")
  }
  # Idempotency: makes the job safe to re-run, and safe to run before a rerun
  # of outstanding jobs lands.
  #
  # Deliberately NOT keyed on BLP_whole. run_blp_whole() returns NULL, by
  # design, whenever tau is constant/near-constant (bug L - GenericML::BLP()
  # cannot identify beta.2), and assigning NULL into a list removes the
  # element - so a successfully patched, degenerate-tau file has no BLP_whole
  # to test for, and would be mistaken for an unpatched one on every later
  # pass. run_independence_test_whole() has no such hole: its tryCatch
  # fallback returns a list even on failure, so independence_cate and
  # independence_po are unconditionally present once add_hte_tests() has run,
  # degenerate tau or not.
  if (!is.null(res[[PATCH_MODEL]]$independence_cate) &&
      !is.null(res[[PATCH_MODEL]]$independence_po)) {
    return("already_patched")
  }
  "to_patch"
}

#' Recompute the three HTE-test fields for one result object
#'
#' Inputs are rebuilt exactly as cate_methods() had them when the arm was
#' fitted: X is the covariate matrix (columns 3..p of `data`), and the nuisances
#' are the whole-sample OOB vectors, which is what run_dr_random_forest() would
#' have been handed.
add_hte_tests <- function(res) {
  X <- as.matrix(res$data[, -c(1:2)])
  Y <- res$data$Y
  W <- res$data$W

  res[[PATCH_MODEL]]$BLP_whole <- run_blp_whole(
    Y, W, res$nuisances_rf$W.hat, res$nuisances_rf$Y0.hat, res[[PATCH_MODEL]]$tau)
  res[[PATCH_MODEL]]$independence_cate <- run_independence_test_whole(
    X, res[[PATCH_MODEL]]$tau)
  # Recomputed rather than copied from causal_forest$independence_po, which is
  # the same test on the same po and the same X. Copying would be faster and
  # give an identical answer; recomputing keeps every value traceable to one
  # function call, and the equivalence check in missing/*/README.md asserts the
  # two agree - which is also how a mis-rebuilt X would be caught.
  res[[PATCH_MODEL]]$independence_po <- run_independence_test_whole(
    X, res$nuisances_rf$po)

  res
}

#' Patch one result file in place
#'
#' The write is atomic: saveRDS to <file>.tmp, then rename over the original.
#' A job killed mid-write leaves an orphan .tmp and an intact result, never a
#' truncated one. An orphan cannot be miscounted either - RES_PATTERN is
#' anchored ("^res_sim_\\d+\\.RDS$"), so "res_sim_1.RDS.tmp" matches nothing.
patch_one_file <- function(path, dry_run = FALSE) {
  t0 <- Sys.time()
  out <- function(status, message = NA_character_, blp_null = NA) {
    data.frame(file = basename(path), status = status, message = message,
               blp_null = blp_null,
               seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 3),
               stringsAsFactors = FALSE)
  }

  res <- tryCatch(readRDS(path),
                  error = function(e) e)
  if (inherits(res, "error")) return(out("error", conditionMessage(res)))

  status <- patch_status_of(res)
  if (status == "already_patched") {
    # Derived from `res`, already in hand - not defaulted to NA. Otherwise a
    # combo's manifest CSV, which is overwritten (not appended) on every call,
    # would lose a real TRUE the moment the job is invoked a second time over
    # an already-patched file, even though nothing changed on disk.
    return(out(status, blp_null = is.null(res[[PATCH_MODEL]]$BLP_whole)))
  }
  if (status != "to_patch") return(out(status))
  if (dry_run) return(out("to_patch"))

  patched <- tryCatch(add_hte_tests(res), error = function(e) e)
  if (inherits(patched, "error")) return(out("error", conditionMessage(patched)))

  tmp <- paste0(path, ".tmp")
  ok <- tryCatch({
    saveRDS(patched, tmp)
    file.rename(tmp, path)
  }, error = function(e) e)

  if (inherits(ok, "error") || !isTRUE(ok)) {
    unlink(tmp)
    msg <- if (inherits(ok, "error")) conditionMessage(ok) else "file.rename returned FALSE"
    return(out("error", msg))
  }
  # Recorded because patch_status_of() can no longer tell from the file alone
  # (that is the point of the fix above) - this is how "how many fits were
  # degenerate?" becomes a manifest query instead of a re-read of 19,800 files.
  out("patched", blp_null = is.null(patched[[PATCH_MODEL]]$BLP_whole))
}

#' Patch a study's results
#'
#' @param study a study_config() object
#' @param combo_idx row(s) of combos(study) to process; NULL means all of them.
#'   The PBS array index is one such row - see missing/*/jobscripts/*_patch.sh.
#' @param dry_run report what would happen; write neither results nor manifest
#' @param verbose one progress line per combination
#' @return invisibly, one row per file processed
patch_hte_tests <- function(study, combo_idx = NULL, dry_run = FALSE,
                            verbose = TRUE) {

  cmb <- combos(study)
  if (is.null(combo_idx)) combo_idx <- seq_len(nrow(cmb))
  stopifnot(all(combo_idx >= 1), all(combo_idx <= nrow(cmb)))

  rows <- bind_rows(lapply(combo_idx, function(k) {
    combo <- cmb[k, , drop = FALSE]
    dir <- combo_dir(study, combo)
    files <- list.files(dir, pattern = RES_PATTERN, full.names = TRUE)
    files <- files[order(run_numbers(files))]

    if (length(files) == 0) {
      if (verbose) cat("combo", k, "-", dir, "- no result files\n")
      return(NULL)
    }

    res_rows <- bind_rows(lapply(files, patch_one_file, dry_run = dry_run))
    res_rows <- bind_cols(combo[rep(1, nrow(res_rows)), , drop = FALSE],
                          tibble(combo_idx = k),
                          res_rows)

    if (verbose) {
      tally <- table(res_rows$status)
      cat("combo", k, "-", nrow(res_rows), "files:",
          paste(names(tally), tally, sep = "=", collapse = " "), "\n")
    }
    res_rows
  }))

  if (is.null(rows) || nrow(rows) == 0) {
    warning("no result files found for the requested combinations")
    return(invisible(rows))
  }

  if (!dry_run) {
    mdir <- patch_manifest_dir(study)
    dir.create(mdir, recursive = TRUE, showWarnings = FALSE)
    # One manifest per combination. A single shared file would be written
    # concurrently by up to 99 array elements and lose rows.
    for (k in unique(rows$combo_idx)) {
      write.csv(rows[rows$combo_idx == k, ],
                file.path(mdir, paste0(k, ".csv")), row.names = FALSE)
    }
  }

  cat("\n", study$name, if (dry_run) " (DRY RUN)" else "", "\n", sep = "")
  print(table(rows$status))
  if (any(rows$status == "error")) {
    warning(sum(rows$status == "error"), " file(s) errored - see the manifest")
  }

  invisible(rows)
}
