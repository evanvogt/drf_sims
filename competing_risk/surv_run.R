##########
# title: run competing-risk simulations 101-500 inside one RStudio session
##########
# The no-queue alternative to `qsub jobscripts/surv_1.sh`. Same grid, same
# results, but run directly inside an interactive session rather than
# waiting on array-job scheduling.
#
# ---- how to run it ----------------------------------------------------------
#   request an RStudio session with 8 cores and 64gb, then
#   source(here::here("competing_risk", "surv_run.R"))
#
# Each grid row runs as its own Rscript subprocess, executing exactly the
# command surv_1.sh runs for one array index:
#
#   Rscript surv_analysis.R <i>
#
# so a row produced here and a row produced by the array are the same
# calculation. surv_analysis.R is not modified, is called with the same
# single argument, and gets the same fresh process a compute node gives it -
# no RNG stream, plan() or OMP_NUM_THREADS state carried over from the row
# before.
#
# ---- workers = 4, not 8 ------------------------------------------------------
# surv_analysis.R hardcodes its own internal `workers <- 2` (the
# plan(multisession, workers = workers) it uses for data generation) - not
# exposed as a CLI arg, so it cannot be lowered from here without editing
# that script. jobscripts/surv_1.sh matches it with ncpus=2 per array task.
# 8 concurrent rows here would therefore ask for 8 x 2 = 16 cores against the
# 8 the session actually has - a 2x oversubscription that measures contention
# rather than running faster. 4 rows x 2 cores each is the true 8-core budget.
# Memory is not the constraint (surv_1.sh asks mem=2gb per task; 4 rows is
# ~8GB against a 64GB session).
#
# ---- default target is runs 101-500, not "every missing row" ----------------
# Unlike cts_val_run.R/me_strategies_run.R, which default to "every row
# without a results file", this script exists specifically to backfill the
# newly widened grid (surv_config.R's run = 1:500, was 1:100) without
# resweeping the original 100 runs per combo, which already have results.
# grid_indices(study, run = 101:500) selects exactly those rows by row
# number - safe, because it does not filter study$grid (R/pipeline.R's
# grid_indices() docs: filtering the grid itself renumbers every row).
#
# Rows whose results file already exists are still skipped (intersected with
# check_failed()), so a session that is interrupted - or hits its walltime -
# is resumed by sourcing this file again.
#
# ---- what it deliberately does not do ---------------------------------------
# It never writes jobscripts/failed_ids.txt and never edits surv_rerun.sh.
# Those belong to the queue path and are surv_check.R's to write; writing
# them from here would leave a resubmission list pointing at rows this
# script has already retried itself. Failed rows are reported and picked up
# by the next source() instead.

library(future)
library(here)

# study$grid, study$res_path, and (via R/pipeline.R) check_failed() and
# grid_indices(). The grid is read from the one place that defines it, so the
# index this script hands surv_analysis.R means what it means everywhere else.
source(here("competing_risk", "surv_config.R"))

# ---- knobs ------------------------------------------------------------------
workers       <- 4      # grid rows in flight at once - 2 cores each, see above
overwrite     <- FALSE  # TRUE re-runs rows that already have a results file
row_timeout   <- 1800   # seconds before a row is killed; 0 disables
poll_interval <- 2      # seconds between checks for a finished row

# EDIT ME to run a different subset - e.g. ids <- grid_indices(study, run = 1)
# for a single run across every combo, or ids <- 1401 for a single row. NULL
# means "runs 101-500, minus whatever already has a results file". Never
# filter study$grid to make a subset: the index is a row number, so filtering
# renumbers every row (R/pipeline.R). Trailing arguments override this, so
# the same file also runs non-interactively: Rscript surv_run.R 1401 1402
ids <- NULL

# ---- what to run --------------------------------------------------------
cli_ids <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)))
if (length(cli_ids) > 0 && !anyNA(cli_ids)) {
  ids <- cli_ids
} else if (is.null(ids)) {
  target <- grid_indices(study, run = 101:500)
  ids <- if (overwrite) target else intersect(target, check_failed(study, write = FALSE))
} else if (!overwrite) {
  # an explicit subset still honours the skip rule, so re-sourcing after an
  # interruption retries only what is actually missing from that subset
  ids <- intersect(ids, check_failed(study, write = FALSE))
}

if (length(ids) && (min(ids) < 1 || max(ids) > nrow(study$grid))) {
  stop("ids must be row numbers of study$grid (1 to ", nrow(study$grid), ")",
       call. = FALSE)
}

# ---- guards -----------------------------------------------------------------
# A row asking for more cores than the session has would measure contention
# rather than run faster - the same check cts_val_run.R/me_strategies_run.R
# make. surv_analysis.R's own workers <- 2 is fixed, so it is baked into
# inner_workers here rather than read from a knob above.
inner_workers <- 2
requested <- workers * inner_workers
available_cores <- as.integer(parallelly::availableCores())
if (requested > available_cores) {
  stop("this session sees ", available_cores, " core(s) but ", workers,
       " rows x ", inner_workers, " core(s) each wants ",
       requested, " - lower `workers`.", call. = FALSE)
}

analysis <- here("competing_risk", "surv_analysis.R")
if (!file.exists(analysis)) {
  stop("cannot find surv_analysis.R at ", analysis, call. = FALSE)
}

# R.home(), not "Rscript" off PATH: rows must run under the same R - and so the
# same library - as the session that launched them, which on the cluster is the
# module/conda environment the array jobs use.
rscript <- file.path(R.home("bin"),
                     if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
if (!file.exists(rscript)) stop("no Rscript at ", rscript, call. = FALSE)

# beside the results tree, outside the repo
log_dir <- file.path(study$res_path, "session_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# ---- the pool ---------------------------------------------------------------
# Each future does nothing but wait on its subprocess, so the 4 worker sessions
# are idle bookkeeping - the 8 cores go to the Rscript children, 2 each.
run_rows <- function(ids) {

  total <- length(ids)
  cat(sprintf("running %d of %d grid rows, %d at a time (%d core(s) available)\n",
              total, nrow(study$grid), workers, available_cores))
  cat(sprintf("logs: %s\n\n", log_dir))

  pb <- utils::txtProgressBar(min = 0, max = total, style = 3)
  on.exit(close(pb), add = TRUE)

  plan(multisession, workers = workers)
  # inside a function this does fire, unlike a top-level loop that has to
  # call plan(sequential) explicitly
  on.exit(plan(sequential), add = TRUE)

  start_row <- function(i) {
    log_file <- file.path(log_dir, paste0("row_", i, ".log"))
    f <- future({
      # captured rather than redirected to a file: on unix, system2() given the
      # same file for stdout and stderr opens it twice and the two streams
      # overwrite each other. suppressWarnings() because a non-zero exit raises
      # one, and the status is handled below rather than as a warning.
      out <- suppressWarnings(
        system2(rscript, c(shQuote(analysis), i),
                stdout = TRUE, stderr = TRUE, timeout = row_timeout))
      writeLines(out, log_file)
      status <- attr(out, "status")
      if (is.null(status)) 0L else as.integer(status)
    }, seed = FALSE)
    list(id = i, started = Sys.time(), f = f, log = log_file)
  }

  queue <- ids
  running <- list()
  done <- 0L
  res_id <- integer(total)
  res_status <- integer(total)
  res_elapsed <- numeric(total)

  while (length(queue) > 0 || length(running) > 0) {

    while (length(running) < workers && length(queue) > 0) {
      running[[length(running) + 1L]] <- start_row(queue[1])
      queue <- queue[-1]
    }

    ready <- which(vapply(running, function(r) resolved(r$f), logical(1)))
    if (length(ready) == 0) {
      Sys.sleep(poll_interval)
      next
    }

    for (r in running[ready]) {
      elapsed <- as.numeric(difftime(Sys.time(), r$started, units = "secs"))
      # a worker that dies takes its future's value with it; that is a failed
      # row, not a reason to abandon the rest of the batch
      status <- tryCatch(value(r$f), error = function(e) {
        message(sprintf("        worker error on row %d: %s", r$id,
                        conditionMessage(e)))
        -1L
      })

      done <- done + 1L
      res_id[done] <- r$id
      res_status[done] <- status
      res_elapsed[done] <- elapsed

      param <- study$grid[r$id, ]
      outcome <- if (status == 0L) {
        sprintf("ok   in %5.1f min", elapsed / 60)
      } else if (status == 124L) {
        sprintf("TIMED OUT after %.0f s - see %s", elapsed, basename(r$log))
      } else {
        sprintf("FAILED status=%d - see %s", status, basename(r$log))
      }

      cat(sprintf("[%4d/%d] row %4d (scenario %g censoring %s run %3d)  %s\n",
                  done, total, r$id, param$scenario, param$censoring,
                  param$run, outcome))
      utils::setTxtProgressBar(pb, done)

      # mean-based, so it is an estimate rather than a promise - but every row
      # here shares the same n and horizon, so the spread is small
      remaining <- total - done
      if (remaining > 0 && done >= workers) {
        eta_min <- mean(res_elapsed[seq_len(done)]) * remaining / workers / 60
        if (done %% workers == 0) {
          cat(sprintf("          ~%.1f h left over %d rows, at the mean so far\n",
                      eta_min / 60, remaining))
        }
      }
    }
    running <- running[-ready]
  }

  data.frame(id = res_id, status = res_status, elapsed = res_elapsed)
}

# ---- go ---------------------------------------------------------------------
if (length(ids) == 0) {
  cat("runs 101-500 already have a results file for every combo - nothing to run.\n")
} else {
  started <- Sys.time()
  outcome <- run_rows(ids)
  wall <- as.numeric(difftime(Sys.time(), started, units = "hours"))

  ok <- sum(outcome$status == 0L)
  timed_out <- sum(outcome$status == 124L)
  failed <- outcome$id[outcome$status != 0L]

  cat(sprintf("\n%d of %d rows ok in %.2f h (%d failed, %d of those timed out)\n",
              ok, nrow(outcome), wall, length(failed), timed_out))
  if (length(failed) > 0) {
    cat("failed rows:", paste(failed, collapse = " "), "\n")
    cat("their logs are in", log_dir, "\n")
  }

  # prints "All simulations complete!" itself when nothing is missing
  missing_after <- intersect(grid_indices(study, run = 101:500),
                              check_failed(study, write = FALSE))
  if (length(missing_after) > 0) {
    cat(sprintf("%d row(s) in runs 101-500 still have no results file - ",
                length(missing_after)),
        "source this file again to retry them\n", sep = "")
  } else {
    cat("every row in runs 101-500 has a results file - nothing missing.\n")
  }
  cat("next: Rscript competing_risk/surv_check.R, then surv_collect.R / surv_metrics.R\n")
}
