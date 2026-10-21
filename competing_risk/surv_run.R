##########
# title: run competing-risk simulations inside one RStudio session
##########
# The no-queue alternative to `qsub jobscripts/surv_1.sh` / `surv_2.sh`: same
# grid, same results, in an interactive session (8 cores, 64gb):
#   source(here::here("competing_risk", "surv_run.R"))
#
# Each row runs as its own `Rscript surv_analysis.R <i>` subprocess, exactly as
# the array does, so no RNG, plan() or thread state carries between rows.
#
# workers = 4: each row uses 2 cores (surv_analysis.R's default workers, as
# surv_1.sh's ncpus=2), so 4 rows fill 8 cores.
#
# By default it targets every row without a results file. The grid is 14,000
# rows, so pick a subset below (e.g. grid_indices(study, rho = 0.5), which
# selects rows without filtering the grid). Re-sourcing resumes.
#
# It never writes failed_ids.txt or edits surv_rerun.sh (study.R check's job),
# so a resubmission list never points at rows this script already retried.

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
# for a single run across every combo, or ids <- 7003 for a single row. NULL
# means "every row without a results file". Never
# filter study$grid to make a subset: the index is a row number, so filtering
# renumbers every row (R/pipeline.R). Trailing arguments override this, so
# the same file also runs non-interactively: Rscript surv_run.R 7003 7004
ids <- NULL

# ---- what to run --------------------------------------------------------
cli_ids <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)))
if (length(cli_ids) > 0 && !anyNA(cli_ids)) {
  ids <- cli_ids
} else if (is.null(ids)) {
  target <- seq_len(nrow(study$grid))
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
# make. surv_analysis.R's default workers = 2 is what every row runs with
# (see the header), so it is baked into inner_workers here rather than read
# from a knob above.
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
  cat("every requested row already has a results file - nothing to run.\n")
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
  missing_after <- intersect(ids, check_failed(study, write = FALSE))
  if (length(missing_after) > 0) {
    cat(sprintf("%d requested row(s) still have no results file - ",
                length(missing_after)),
        "source this file again to retry them\n", sep = "")
  } else {
    cat("every requested row has a results file - nothing missing.\n")
  }
  cat("next: Rscript study.R check competing_risk/surv_config.R, then collect / surv_metrics.R\n")
}
