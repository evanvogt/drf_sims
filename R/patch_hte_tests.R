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

#' How often a combination's manifest is flushed to disk mid-run
#'
#' The manifest used to be written only once every file in the combination had
#' been processed, which meant an array element killed part-way through left no
#' manifest at all - and so was invisible to check_all.R, which counts manifest
#' rows. That is exactly how ten of missing/binary's 99 elements went missing
#' without a trace; see missing/binary/bin_miss_patch_check.R. Flushing
#' periodically costs one small write.csv per 25 files and turns a killed
#' element into a short manifest that says how far it got.
MANIFEST_FLUSH_EVERY <- 25L

#' Write one combination's manifest CSV
#'
#' One file per combination, never a shared one: this is written concurrently by
#' up to 99 array elements, and a single shared file would lose rows.
#'
#' Split out so the periodic flush and the final write cannot drift in shape -
#' a partial manifest must have exactly the columns of a complete one, or
#' check_all.R's bind_rows() over the directory fails on the mix.
write_manifest <- function(mdir, k, combo, res_rows) {
  rows <- bind_cols(combo[rep(1, nrow(res_rows)), , drop = FALSE],
                    tibble(combo_idx = k),
                    res_rows)
  write.csv(rows, file.path(mdir, paste0(k, ".csv")), row.names = FALSE)
  rows
}

#' Patch a study's results
#'
#' @param study a study_config() object
#' @param combo_idx row(s) of combos(study) to process; NULL means all of them.
#'   The PBS array index is one such row - see missing/*/jobscripts/*_patch.sh.
#' @param dry_run report what would happen; write neither results nor manifest
#' @param verbose one progress line per FILE, plus a tally per combination. The
#'   per-file lines go to stdout on purpose: the manifest is only flushed every
#'   MANIFEST_FLUSH_EVERY files, so between flushes the job's .o log is the only
#'   record of where it had got to when it died.
#' @return invisibly, one row per file processed
patch_hte_tests <- function(study, combo_idx = NULL, dry_run = FALSE,
                            verbose = TRUE) {

  cmb <- combos(study)
  if (is.null(combo_idx)) combo_idx <- seq_len(nrow(cmb))
  stopifnot(all(combo_idx >= 1), all(combo_idx <= nrow(cmb)))

  # Created up front rather than after the loop, because the flush inside it
  # needs somewhere to write.
  mdir <- patch_manifest_dir(study)
  if (!dry_run) dir.create(mdir, recursive = TRUE, showWarnings = FALSE)

  rows <- bind_rows(lapply(combo_idx, function(k) {
    combo <- cmb[k, , drop = FALSE]
    dir <- combo_dir(study, combo)
    files <- list.files(dir, pattern = RES_PATTERN, full.names = TRUE)
    files <- files[order(run_numbers(files))]

    if (length(files) == 0) {
      if (verbose) cat("combo", k, "-", dir, "- no result files\n")
      return(NULL)
    }

    acc <- vector("list", length(files))
    for (fi in seq_along(files)) {
      acc[[fi]] <- patch_one_file(files[fi], dry_run = dry_run)
      if (verbose) {
        cat("  combo", k, "file", basename(files[fi]), "-",
            acc[[fi]]$status, "\n")
      }
      if (!dry_run && fi %% MANIFEST_FLUSH_EVERY == 0 && fi < length(files)) {
        write_manifest(mdir, k, combo, bind_rows(acc[seq_len(fi)]))
      }
    }

    # The complete write. Also the only write when the combination is shorter
    # than one flush interval.
    res_rows <- bind_rows(acc)
    res_rows <- if (dry_run) {
      bind_cols(combo[rep(1, nrow(res_rows)), , drop = FALSE],
                tibble(combo_idx = k), res_rows)
    } else {
      write_manifest(mdir, k, combo, res_rows)
    }

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

  cat("\n", study$name, if (dry_run) " (DRY RUN)" else "", "\n", sep = "")
  print(table(rows$status))
  if (any(rows$status == "error")) {
    warning(sum(rows$status == "error"), " file(s) errored - see the manifest")
  }

  invisible(rows)
}

##### auditing the repair #####################################################
# check_all.R reports a study's patch progress by counting manifest rows, which
# is fast and, when an array element dies, wrong in a way that hides the
# failure: no manifest is written until MANIFEST_FLUSH_EVERY files are done, so
# a killed element looks identical to one that never ran. missing/binary sat at
# 7,800/8,800 for exactly this reason - ten of its 99 elements left no manifest.
#
# Reading the job's stderr is the obvious next move and was not available: the
# HPC returns no .e files and the .o files from that run are gone. Everything
# below therefore reconstructs the diagnosis from the result files, which do
# survive. Three traces make that possible:
#
#   1. patch_status_of() is authoritative about whether a file was patched;
#   2. mtimes say WHEN, because patch_one_file() rewrites the file - so the
#      first and last patched file bracket the element's runtime, and files it
#      never reached still carry their much older timestamp from the simulation;
#   3. an orphan res_sim_<n>.RDS.tmp names the file it was writing when it died.

#' Orphan half-writes left by a killed patch job
#'
#' The mirror of RES_PATTERN, and anchored for the same reason: patch_one_file()
#' writes <file>.tmp and renames, so "res_sim_1.RDS.tmp" is a job that died
#' between the two. RES_PATTERN deliberately does not match it (it would inflate
#' every file count in the repo); this exists to find it on purpose.
TMP_PATTERN <- "^res_sim_\\d+\\.RDS\\.tmp$"

#' Which files belong to the most recent burst of writes
#'
#' Used only when the caller cannot afford to open every file. The patch
#' rewrites files in run order over a few minutes, so its files cluster tightly
#' in mtime while the ones it never reached sit hours or days earlier, at
#' whenever the simulation wrote them. Split on the last gap of at least
#' `min_gap_secs` and take everything after it.
#'
#' Deliberately fallible, which is why patch_forensics() labels its output
#' `basis = "mtime"`: a combination NOTHING patched has no gap either, and is
#' indistinguishable here from one that was patched end to end. Only the
#' `deep = TRUE` path can tell those apart, and it is the default.
#'
#' @return logical, aligned to `mtime`; TRUE = in the most recent burst
mtime_window <- function(mtime, min_gap_secs = 3600) {
  if (length(mtime) < 2 || anyNA(mtime)) return(rep(TRUE, length(mtime)))
  ord <- order(mtime)
  m <- mtime[ord]
  gaps <- as.numeric(difftime(m[-1], m[-length(m)], units = "secs"))
  big <- which(gaps >= min_gap_secs)
  if (!length(big)) return(rep(TRUE, length(mtime)))
  cut <- max(big)
  out <- logical(length(mtime))
  out[ord] <- c(rep(FALSE, cut), rep(TRUE, length(m) - cut))
  out
}

#' Manifest-level status of every parameter combination
#'
#' Fast - one list.files() and one read.csv() per combination, no result file is
#' opened. `done` uses the same definition as check_all.R's count_patched()
#' (patched + already_patched) so the totals here reconcile with the campaign
#' table rather than offering a second, subtly different number.
#'
#' multiple_imputation is excluded from the denominator exactly as count_patched()
#' excludes it: those runs keep no nuisances, the patch refuses them by design,
#' and counting them would peg a fully repaired study at 88.9%. See the
#' multiple-imputation note in missing/README.md.
patch_manifest_status <- function(study) {
  cmb <- combos(study)
  mdir <- patch_manifest_dir(study)

  bind_rows(lapply(seq_len(nrow(cmb)), function(k) {
    combo <- cmb[k, , drop = FALSE]
    n_files <- length(list.files(combo_dir(study, combo), pattern = RES_PATTERN))
    mf <- file.path(mdir, paste0(k, ".csv"))

    rows <- if (file.exists(mf)) {
      tryCatch(read.csv(mf, stringsAsFactors = FALSE), error = function(e) NULL)
    } else NULL

    tally <- function(s) if (is.null(rows)) 0L else sum(rows$status == s)
    done <- if (is.null(rows)) 0L else {
      sum(rows$status %in% c("patched", "already_patched"))
    }

    is_mi <- "method" %in% names(combo) &&
      identical(as.character(combo$method), "multiple_imputation")
    expected_done <- if (is_mi) 0L else n_files

    errs <- if (is.null(rows)) character(0) else {
      unique(rows$message[rows$status == "error"])
    }
    errs <- errs[!is.na(errs)]

    bind_cols(
      tibble(combo_idx = k),
      as_tibble(combo),
      tibble(
        is_mi            = is_mi,
        n_files          = n_files,
        manifest         = file.exists(mf),
        manifest_rows    = if (is.null(rows)) NA_integer_ else nrow(rows),
        done             = done,
        expected_done    = expected_done,
        shortfall        = max(expected_done - done, 0L),
        n_skipped_no_arm = tally("skipped_no_arm"),
        n_error          = tally("error"),
        n_to_patch       = tally("to_patch"),
        error_message    = if (length(errs)) {
          paste(errs, collapse = " | ")
        } else NA_character_))
  }))
}

#' On-disk forensics for a set of parameter combinations
#'
#' @param combo_idx rows of combos(study) to examine
#' @param deep open every file and apply patch_status_of(). Authoritative, and
#'   about 100 readRDS calls per combination - fine for the handful that look
#'   wrong, too slow for all 99. FALSE falls back to mtime_window(), which is
#'   cheap enough to run over the whole study and is what supplies the
#'   comparison timings patch_failure_signature() needs.
patch_forensics <- function(study, combo_idx, deep = TRUE) {
  cmb <- combos(study)

  bind_rows(lapply(combo_idx, function(k) {
    combo <- cmb[k, , drop = FALSE]
    dir <- combo_dir(study, combo)
    files <- list.files(dir, pattern = RES_PATTERN, full.names = TRUE)
    runs <- run_numbers(files)
    ord <- order(runs)
    files <- files[ord]
    runs <- runs[ord]
    orphans <- list.files(dir, pattern = TMP_PATTERN)

    out <- tibble(
      combo_idx           = k,
      basis               = NA_character_,
      n_orphan_tmp        = length(orphans),
      orphan_tmp          = if (length(orphans)) {
        paste(orphans, collapse = ",")
      } else NA_character_,
      n_patched           = NA_integer_,
      n_unpatched         = NA_integer_,
      n_skipped_mi_files  = NA_integer_,
      n_no_arm_files      = NA_integer_,
      n_unreadable        = NA_integer_,
      first_unpatched_run = NA_integer_,
      last_patched_run    = NA_integer_,
      patch_window_start  = as.POSIXct(NA),
      died_at             = as.POSIXct(NA),
      span_minutes        = NA_real_)

    if (!length(files)) return(out)

    mtime <- file.info(files)$mtime

    if (deep) {
      st <- vapply(files, function(f) {
        res <- tryCatch(readRDS(f), error = function(e) e)
        if (inherits(res, "error")) "unreadable" else patch_status_of(res)
      }, character(1), USE.NAMES = FALSE)
      out$basis <- "files"
      patched <- st == "already_patched"
      out$n_skipped_mi_files <- sum(st == "skipped_mi")
      out$n_no_arm_files     <- sum(st == "skipped_no_arm")
      out$n_unreadable       <- sum(st == "unreadable")
      unpatched <- st == "to_patch"
    } else {
      out$basis <- "mtime"
      patched <- mtime_window(mtime)
      unpatched <- !patched
    }

    out$n_patched   <- sum(patched)
    out$n_unpatched <- sum(unpatched)
    if (any(unpatched)) out$first_unpatched_run <- min(runs[unpatched])
    if (any(patched)) {
      out$last_patched_run   <- max(runs[patched])
      out$patch_window_start <- min(mtime[patched])
      out$died_at            <- max(mtime[patched])
      out$span_minutes <- round(
        as.numeric(difftime(out$died_at, out$patch_window_start,
                            units = "mins")), 2)
    }
    out
  }))
}

#' Walltime of a patch jobscript, in minutes
#'
#' Reuses pipeline.R's read_script()/parse_walltime() rather than re-reading the
#' PBS directives a third way. NULL when the script is absent or unparseable -
#' the signature below then simply drops the walltime rule.
patch_walltime_minutes <- function(script_path) {
  if (is.null(script_path) || !file.exists(script_path)) return(NULL)
  line <- grep("^#PBS -l walltime=", read_script(script_path)$lines, value = TRUE)
  if (!length(line)) return(NULL)
  secs <- parse_walltime(line[1])
  if (is.null(secs)) NULL else secs / 60
}

#' Why the failed elements failed, inferred from timings alone
#'
#' This is the part that stands in for the missing logs. Every element ran the
#' same code over the same number of files, so the combinations that finished
#' calibrate what "normal" costs; the ones that did not are then read against
#' that baseline and against the job's own walltime.
#'
#' @param fx forensics for EVERY combination (deep or not - only the timings and
#'   orphan counts are used)
#' @param suspects combo indices that failed
#' @param walltime_minutes from patch_walltime_minutes(), or NULL
#' @return list(verdict, detail) - detail is character lines for the report
patch_failure_signature <- function(fx, suspects, walltime_minutes = NULL) {
  ok  <- fx[!fx$combo_idx %in% suspects & !is.na(fx$span_minutes), ]
  bad <- fx[fx$combo_idx %in% suspects, ]
  timed <- bad[!is.na(bad$span_minutes), ]

  detail <- character(0)
  typical <- if (nrow(ok)) median(ok$span_minutes) else NA_real_
  if (!is.na(typical)) {
    detail <- c(detail, sprintf(
      "Combinations that completed: %d, median span %.1f min (max %.1f).",
      nrow(ok), typical, max(ok$span_minutes)))
  }
  if (!nrow(timed)) {
    return(list(verdict = "no timing evidence - the failed combinations show no patched files at all",
                detail = detail))
  }

  detail <- c(detail, sprintf(
    "Failed combinations with partial work: %d, median span %.1f min (max %.1f).",
    nrow(timed), median(timed$span_minutes), max(timed$span_minutes)))
  if (any(bad$n_orphan_tmp > 0)) {
    detail <- c(detail, sprintf(
      "Orphan .tmp files in %d of them - killed mid-write, so the kill was external.",
      sum(bad$n_orphan_tmp > 0)))
  }
  # How tightly the stops cluster in wall-clock time. A common instant points at
  # something that hit the whole job at once; a common DURATION points at each
  # element hitting its own limit.
  spread <- if (nrow(timed) > 1) {
    as.numeric(difftime(max(timed$died_at), min(timed$died_at), units = "mins"))
  } else NA_real_
  if (!is.na(spread)) {
    detail <- c(detail, sprintf(
      "They stopped within a %.1f-minute window of each other (%s to %s).",
      spread, format(min(timed$died_at), "%Y-%m-%d %H:%M"),
      format(max(timed$died_at), "%H:%M")))
  }
  if (!is.null(walltime_minutes)) {
    detail <- c(detail, sprintf("Job walltime: %.0f min.", walltime_minutes))
  }

  verdict <- if (!is.null(walltime_minutes) &&
                 median(timed$span_minutes) >= 0.85 * walltime_minutes) {
    paste0("killed on walltime - the failed elements ran to ~",
           round(median(timed$span_minutes)), " min against a ",
           round(walltime_minutes), "-minute limit. Raise the walltime and ",
           "throttle the array so the elements are not all hammering the ",
           "filesystem at once.")
  } else if (!is.na(spread) && spread <= 5 && nrow(timed) > 1) {
    paste0("stopped together within ", round(spread, 1), " minutes, well short ",
           "of the walltime - a node, filesystem or scheduler event took them ",
           "out rather than each hitting its own limit.")
  } else if (any(bad$n_orphan_tmp > 0)) {
    "killed mid-write at scattered times - external kill, cause not identifiable from timings alone."
  } else if (!is.na(typical) && median(timed$span_minutes) < typical) {
    paste0("stopped early, at scattered times, with no half-written file - the ",
           "signature of an R-level error rather than a kill. The re-run will ",
           "surface it now that stderr is captured.")
  } else {
    "inconclusive from timings - read the re-run's log."
  }

  list(verdict = verdict, detail = detail)
}

#' Which parameter combinations the patch failed on, and why
#'
#' The counterpart to pipeline.R's check_failed(), for the repair rather than
#' the simulation, and it follows the same loop: report, write a list of indices
#' to resubmit, and leave the qsub to the user. It writes its report next to the
#' study's scripts rather than into the results tree, because results never
#' leave the HPC - a committed CSV/MD is how the answer reaches a laptop.
#'
#' @param deep open every result file in a suspect combination. This is the
#'   check that separates "the element died before writing its manifest but the
#'   files really are patched" from "the files are genuinely unpatched" - the
#'   two need different responses, and the manifest cannot tell them apart.
#' @param write write failed_patch_ids.txt and the report
#' @param report_dir where the CSV/MD go; defaults to the study's script folder
#' @param failed_file where the resubmission list goes
#' @param patch_script the jobscript, read only for its walltime
#' @return invisibly, the per-combination table
check_patch_failed <- function(study, deep = TRUE, write = TRUE,
                               report_dir = NULL, failed_file = NULL,
                               patch_script = NULL, verbose = TRUE) {

  if (is.null(study$failed_file) && (is.null(report_dir) || is.null(failed_file))) {
    stop("study has no failed_file - pass report_dir and failed_file explicitly")
  }
  job_dir <- dirname(study$failed_file)
  if (is.null(report_dir))   report_dir   <- dirname(job_dir)
  if (is.null(failed_file))  failed_file  <- file.path(job_dir, "failed_patch_ids.txt")
  if (is.null(patch_script)) {
    patch_script <- file.path(job_dir, paste0(study$prefix, "_patch.sh"))
  }

  if (verbose) cat("reading manifests...\n")
  st <- patch_manifest_status(study)

  # Cheap pass over everything: the completed combinations are the baseline the
  # failures are judged against, so they have to be measured too.
  if (verbose) cat("timing all", nrow(st), "combinations from file mtimes...\n")
  fx <- patch_forensics(study, st$combo_idx, deep = FALSE)

  # A combination is suspect if it owes patched files, or its manifest is
  # missing/short/errored. MI combinations owe nothing but are still flagged
  # when their manifest is absent, so the record is complete.
  suspect <- st$shortfall > 0 | st$n_error > 0 |
    (!st$manifest & st$n_files > 0) |
    (!is.na(st$manifest_rows) & st$manifest_rows < st$n_files)

  if (deep && any(suspect)) {
    idx <- st$combo_idx[suspect]
    if (verbose) {
      cat("opening result files for", length(idx), "suspect combination(s) -",
          sum(st$n_files[suspect]), "files...\n")
    }
    deep_fx <- patch_forensics(study, idx, deep = TRUE)
    fx <- bind_rows(deep_fx, fx[!fx$combo_idx %in% idx, ])
    fx <- fx[order(fx$combo_idx), ]
  }

  res <- left_join(st, fx, by = "combo_idx")

  res$diagnosis <- vapply(seq_len(nrow(res)), function(i) {
    r <- res[i, ]
    if (r$n_files == 0)          return("no_results")
    if (r$n_error > 0)           return("manifest_errors")

    # Claims about the DATA rather than about a job that died, so they are
    # settled first and from whichever source saw them. Reading them off the
    # manifest alone would be wrong in exactly the case this script exists for:
    # the manifest is missing, and the shortfall it would have explained gets
    # blamed on a kill and pointlessly resubmitted.
    if (r$n_skipped_no_arm > 0 ||
        (!is.na(r$n_no_arm_files) && r$n_no_arm_files > 0)) {
      return("skipped_no_arm")
    }
    # A non-MI combination whose files kept no nuisances can never be patched -
    # add_hte_tests() has nothing to compute a BLP from. Distinct from
    # mi_by_design because there it is expected; here it is a finding.
    if (!r$is_mi && !is.na(r$n_skipped_mi_files) && r$n_skipped_mi_files > 0) {
      return("no_nuisances")
    }

    if (r$is_mi) return(if (r$manifest) "mi_by_design" else "mi_no_manifest")
    if (!r$manifest) {
      if (is.na(r$n_patched))                       return("no_manifest_unknown")
      if (r$n_unpatched == 0 && r$n_patched > 0)    return("no_manifest_files_patched")
      if (r$n_patched == 0 && r$n_unpatched > 0)    return("no_manifest_files_unpatched")
      return("no_manifest_partial")
    }
    if (!is.na(r$manifest_rows) && r$manifest_rows < r$n_files) return("manifest_short")
    if (r$shortfall > 0)         return("incomplete")
    "ok"
  }, character(1))

  # skipped_no_arm and no_nuisances are excluded on purpose: both are claims
  # about the data, not jobs that died, and re-running the patch will not change
  # either. They are still reported - a permanent shortfall that nobody can
  # close is worth knowing about, and is not the same as an unsubmitted job.
  no_rerun <- c("ok", "mi_by_design", "no_results", "skipped_no_arm",
                "no_nuisances")
  rerun <- sort(res$combo_idx[!res$diagnosis %in% no_rerun])

  sig <- patch_failure_signature(fx, rerun, patch_walltime_minutes(patch_script))

  # ---- console ----
  cat("\n=== ", study$name, " patch audit ===\n", sep = "")
  print(table(res$diagnosis))
  cat("\npatchable:", sum(res$expected_done),
      " patched:", sum(res$done),
      " shortfall:", sum(res$shortfall), "\n")

  if (length(rerun)) {
    # Columns that take one value across the whole study (n, type, prop in the
    # missing-data grids) say nothing here and are what pushes the table over
    # the console width into an unreadable three-block wrap.
    varying <- study$path_cols[vapply(study$path_cols, function(cl) {
      length(unique(res[[cl]])) > 1
    }, logical(1))]
    cols <- c("combo_idx", varying, "n_files", "done", "diagnosis",
              "n_patched", "first_unpatched_run", "n_orphan_tmp",
              "span_minutes")
    cat("\n")
    old_width <- options(width = 200)
    on.exit(options(old_width), add = TRUE)
    print(as.data.frame(res[res$combo_idx %in% rerun, cols]), row.names = FALSE)
    cat("\n--- why ---\n")
    for (l in sig$detail) cat(" ", l, "\n")
    cat("\nVERDICT: ", sig$verdict, "\n", sep = "")
    errs <- unique(res$error_message[!is.na(res$error_message)])
    if (length(errs)) {
      cat("\nerror messages recorded in the manifests:\n")
      for (e in errs) cat("  ", e, "\n")
    }
  } else {
    cat("\nNothing to re-run - every patchable combination is accounted for.\n")
  }

  if (write) {
    dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
    csv_path <- file.path(report_dir, paste0(study$prefix, "_patch_check.csv"))
    md_path  <- file.path(report_dir, paste0(study$prefix, "_patch_check.md"))
    write.csv(as.data.frame(res), csv_path, row.names = FALSE)
    writeLines(patch_check_md(study, res, rerun, sig), md_path)

    # -J 1-0 is not a legal array, so an empty list means leaving the file
    # alone rather than writing something unsubmittable - same reasoning as
    # check_failed()'s "nothing to submit" branch.
    if (length(rerun)) {
      cat(rerun, file = failed_file, sep = "\n")
      cat("\nwrote", length(rerun), "index/indices to", failed_file, "\n")
      cat("set -J 1-", length(rerun), " in ", study$prefix,
          "_patch_rerun.sh before submitting\n", sep = "")
    }
    cat("wrote", csv_path, "and", md_path, "\n")
  }

  invisible(res)
}

#' The report check_patch_failed() leaves behind for `git pull` to carry home
patch_check_md <- function(study, res, rerun, sig) {
  cols <- c("combo_idx", study$path_cols, "n_files", "manifest_rows", "done",
            "diagnosis", "n_patched", "first_unpatched_run", "n_orphan_tmp",
            "died_at", "span_minutes", "error_message")
  cols <- intersect(cols, names(res))
  tbl <- as.data.frame(res[res$combo_idx %in% rerun, cols])
  tbl[] <- lapply(tbl, function(x) ifelse(is.na(x), "", as.character(x)))

  md_rows <- if (nrow(tbl)) {
    c(paste0("| ", paste(names(tbl), collapse = " | "), " |"),
      paste0("|", paste(rep("---", ncol(tbl)), collapse = "|"), "|"),
      apply(tbl, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
  } else "Every patchable combination is accounted for."

  c(paste0("# ", study$name, " - HTE patch audit"),
    "",
    paste0("Last updated: ", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
    "",
    sprintf("patchable %d | patched %d | shortfall %d | combinations to re-run %d",
            sum(res$expected_done), sum(res$done), sum(res$shortfall),
            length(rerun)),
    "",
    "## Combinations needing a re-run",
    "",
    md_rows,
    "",
    "## Why they failed",
    "",
    paste0("**", sig$verdict, "**"),
    "",
    paste0("- ", sig$detail),
    "",
    "## Diagnoses",
    "",
    paste0("- `", names(table(res$diagnosis)), "`: ",
           as.integer(table(res$diagnosis))),
    "",
    "Generated by `<prefix>_patch_check.R` -> `check_patch_failed()` in",
    "`R/patch_hte_tests.R`. Re-run the listed indices with",
    paste0("`qsub jobscripts/", study$prefix, "_patch_rerun.sh`; the patch is"),
    "idempotent, so re-running a combination that is already correct is safe.")
}
