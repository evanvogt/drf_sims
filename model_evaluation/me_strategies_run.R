##########
# title: run the me_strategies nuisance-arm pass inside one RStudio session
##########
# The no-queue alternative to `qsub jobscripts/me_strategies.sh`. Same 360
# grid rows and the same results, but 6 at a time inside an interactive
# session rather than 360 array jobs waiting to be scheduled.
#
# ---- how to run it ----------------------------------------------------------
#   request an RStudio session with 8 cores and 64gb, then
#   source(here::here("model_evaluation", "me_strategies_run.R"))
#
# Each grid row runs as its own Rscript subprocess, executing exactly the
# command me_strategies.sh runs for one array index:
#
#   Rscript me_strategies.R <i> 1
#
# so a row produced here and a row produced by the array are the same
# calculation. me_strategies.R is not modified, is called with the same
# arguments, and gets the same fresh process a compute node gives it - no RNG
# kind, plan() or OMP_NUM_THREADS state carried over from the row before.
#
# ---- workers = 6, not 8 ------------------------------------------------------
# me_strategies.R hardcodes each row's H2O JVM heap to max_mem_size = "10G"
# (me_strategies.R:49) - not exposed as a CLI arg, so it cannot be lowered
# from here without editing that script. That cap is also why the PBS array
# job (me_strategies.sh) throttles to 4 concurrent tasks despite requesting
# only 2 cores each. 8 concurrent rows on a 64GB session could ask for up to
# 8 x 10G = 80GB of heap alone, before the R/XGBoost processes or the OS; 6
# x 10G = 60GB leaves headroom. Each row still gets a single core
# (n_cores = 1) - only the number running at once is capped below 8.
#
# ---- the 2 permanently-missing rows ------------------------------------------
# 2 of the 360 grid rows failed repeatedly in the main study and have no
# me_analysis.R source file (see me_check.R's header). me_strategies.R exits
# 0 for them without writing a result - there is nothing for it to read. A
# plain "which study_strat rows are missing a results file" check would
# therefore report those 2 rows missing forever, so this script uses
# reachable_missing() below (the same logic me_check.R uses) to exclude rows
# that are missing only because their source run is itself missing. A
# complete pass reports 358/358, not 358/360.
#
# ---- what it deliberately does not do ---------------------------------------
# It never writes jobscripts/failed_ids_strat.txt and never edits
# me_strat_rerun.sh. Those belong to the queue path and are me_check.R's to
# write; writing them from here would leave a resubmission list pointing at
# rows this script has already retried itself. Failed rows are reported and
# picked up by the next source() instead.

library(future)
library(here)

# study (for the reachable-missing check) and study_strat (the target of
# this pass) come from the one place that defines the grid, so the index
# this script hands me_strategies.R means what it means everywhere else.
source(here("model_evaluation", "me_config.R"))

# ---- knobs ------------------------------------------------------------------
workers       <- 6      # grid rows in flight at once - one core each. See the
                         # H2O memory note above for why this is 6, not 8.
n_cores       <- 1      # -> me_strategies.R's `n_cores` arg (XGBoost nthread
                         # / H2O nthreads) - the single core per run asked for
overwrite     <- FALSE  # TRUE re-runs rows that already have a results file
row_timeout   <- 1800   # seconds before a row is killed; 0 disables
poll_interval <- 2      # seconds between checks for a finished row

# EDIT ME to run a subset - e.g. ids <- grid_indices(study_strat, scenario = 1)
# for one scenario, or ids <- 1 for a single row. NULL means "every
# study_strat row still reachably missing" (see reachable_missing() below).
# Never filter study_strat$grid to make a subset: the index is a row number,
# so filtering renumbers every row (R/pipeline.R). Trailing arguments
# override this, so the same file also runs non-interactively:
# Rscript me_strategies_run.R 1 2 3
ids <- NULL

# ---- reachable-missing helper -----------------------------------------------
# A study_strat row only counts as missing if its (scenario, n, run) key is
# not ALSO missing from `study` (i.e. it has a me_analysis.R source file to
# read). See the header note above.
reachable_missing <- function() {
  missing_idx <- check_failed(study_strat, write = FALSE)
  src_missing <- check_failed(study, write = FALSE)
  key <- function(st, idx) {
    do.call(paste, c(st$grid[idx, c(st$path_cols, "run")], sep = "/"))
  }
  missing_idx[!(key(study_strat, missing_idx) %in% key(study, src_missing))]
}

# ---- what to run --------------------------------------------------------
cli_ids <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)))
if (length(cli_ids) > 0 && !anyNA(cli_ids)) {
  ids <- cli_ids
} else if (is.null(ids)) {
  ids <- if (overwrite) seq_len(nrow(study_strat$grid)) else reachable_missing()
} else if (!overwrite) {
  # an explicit subset still honours the skip rule, so re-sourcing after an
  # interruption retries only what is actually missing from that subset
  ids <- intersect(ids, reachable_missing())
}

if (length(ids) && (min(ids) < 1 || max(ids) > nrow(study_strat$grid))) {
  stop("ids must be row numbers of study_strat$grid (1 to ",
       nrow(study_strat$grid), ")", call. = FALSE)
}

# ---- guards -----------------------------------------------------------------
# A row asking for more cores than the session has would measure contention
# rather than run faster - the same check cts_val_run.R makes.
requested <- workers * n_cores
available_cores <- as.integer(parallelly::availableCores())
if (requested > available_cores) {
  stop("this session sees ", available_cores, " core(s) but ", workers,
       " rows x ", n_cores, " core(s) each wants ",
       requested, " - lower `workers`.", call. = FALSE)
}

analysis <- here("model_evaluation", "me_strategies.R")
if (!file.exists(analysis)) {
  stop("cannot find me_strategies.R at ", analysis, call. = FALSE)
}

# R.home(), not "Rscript" off PATH: rows must run under the same R - and so the
# same library - as the session that launched them, which on the cluster is the
# module/conda environment the array jobs use.
rscript <- file.path(R.home("bin"),
                     if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
if (!file.exists(rscript)) stop("no Rscript at ", rscript, call. = FALSE)

# beside the results tree, outside the repo
log_dir <- file.path(study_strat$res_path, "session_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# ---- the pool ---------------------------------------------------------------
# Each future does nothing but wait on its subprocess, so the 6 worker sessions
# are idle bookkeeping - the cores go to the Rscript children, one row each.
run_rows <- function(ids) {

  total <- length(ids)
  cat(sprintf("running %d of %d grid rows, %d at a time (%d core(s) available)\n",
              total, nrow(study_strat$grid), workers, available_cores))
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
      # the child inherits this. Left unset, each of the concurrent rows
      # would let its BLAS/OpenMP spawn a thread per visible core - what
      # n_cores = 1 prevents on the cluster.
      Sys.setenv(OMP_NUM_THREADS = n_cores)
      # captured rather than redirected to a file: on unix, system2() given the
      # same file for stdout and stderr opens it twice and the two streams
      # overwrite each other. suppressWarnings() because a non-zero exit raises
      # one, and the status is handled below rather than as a warning.
      out <- suppressWarnings(
        system2(rscript, c(shQuote(analysis), i, n_cores),
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
      # row, not a reason to abandon the other 359
      status <- tryCatch(value(r$f), error = function(e) {
        message(sprintf("        worker error on row %d: %s", r$id,
                        conditionMessage(e)))
        -1L
      })

      done <- done + 1L
      res_id[done] <- r$id
      res_status[done] <- status
      res_elapsed[done] <- elapsed

      param <- study_strat$grid[r$id, ]
      outcome <- if (status == 0L) {
        sprintf("ok   in %5.1f min", elapsed / 60)
      } else if (status == 124L) {
        sprintf("TIMED OUT after %.0f s - see %s", elapsed, basename(r$log))
      } else {
        sprintf("FAILED status=%d - see %s", status, basename(r$log))
      }

      cat(sprintf("[%4d/%d] row %4d (scenario %g n %g run %3d)  %s\n",
                  done, total, r$id, param$scenario, param$n, param$run, outcome))
      utils::setTxtProgressBar(pb, done)

      # mean-based, so it is an estimate rather than a promise - but every row
      # here shares the same shape of work, so the spread is small
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
  cat("every reachable row already has a results file - nothing to run.\n")
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

  missing_after <- reachable_missing()
  if (length(missing_after) > 0) {
    cat(sprintf("%d row(s) still have no results file - source this file again ",
                length(missing_after)),
        "to retry them\n", sep = "")
  } else {
    cat("every reachable row has a results file - nothing missing.\n")
  }
  cat("next: Rscript model_evaluation/me_check.R strategies, then me_collect.R / me_metrics.R\n")
}
