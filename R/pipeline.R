##########
# title: shared simulation-pipeline plumbing
##########
# One copy of get_results() and check_failed(), which existed in eight near-
# identical versions differing only in a parameter grid and a results path.
#
# The point is not the line count. Each study defined its grid FOUR times - in
# *_analysis.R, *_check.R, *_collect.R and implicitly in *_metrics.R - and the
# copies had drifted:
#
#   * missing/continuous/cts_miss_collect.R looked for mechanism "AUX"/"AUX-Y"
#     while the DGM, analysis and check scripts all used "MNAR"/"MNAR-Y", so two
#     of the three mechanisms were silently collected as empty.  (bug B)
#   * cts_miss_metrics.R built its relative-efficiency reference from
#     method == "complete_data", which the collect grid never included, so
#     rel_efficiency was NA for every row.  (bug C)
#   * cts_miss_analysis.R filtered the grid AFTER expand.grid, which renumbered
#     the array index, so index i meant different parameters in the analysis
#     script than in the check and rerun scripts.  (bug D)
#
# A study now declares its grid once, in <study>/<prefix>_config.R, and every
# script reads it from there. The array index is always a row of `grid`.

require(dplyr)
require(purrr)
require(tibble)

#' Declare a study
#'
#' @param name human-readable id, used in messages
#' @param prefix the file-name stem this study's scripts share - "bin" for
#'   bin_config.R, bin_analysis.R, jobscripts/bin_1.sh, jobscripts/bin_rerun.sh.
#'   update_rerun_script() uses it to find the last two. It is declared rather
#'   than inferred because confidence_intervals/optimal_sf/jobscripts holds two
#'   studies, so globbing that directory for *_rerun.sh matches both.
#' @param res_path directory holding this study's results
#' @param grid THE parameter grid. One row per array job, in array-index order,
#'   and it must include a `run` column. Never filter this after the fact - the
#'   PBS array index is a row number, so filtering renumbers every job.
#' @param path_cols the grid columns that form the results subdirectory path,
#'   in order.
#' @param path_prefix named character vector of directory-name prefixes, e.g.
#'   c(scenario = "scenario_", censoring = "censor_"). Columns not named here
#'   are rendered as-is. Defaults to prefixing `scenario`, which every study does.
#' @param n_sims runs per parameter combination
#' @param failed_file where check_failed() writes the resubmission list
study_config <- function(name, prefix, res_path, grid, path_cols, n_sims = 100,
                         failed_file = NULL,
                         path_prefix = c(scenario = "scenario_")) {
  stopifnot("run" %in% names(grid))
  stopifnot(all(path_cols %in% names(grid)))
  stopifnot(is.character(prefix), length(prefix) == 1, nzchar(prefix))
  list(name = name, prefix = prefix, res_path = res_path, grid = grid,
       path_cols = path_cols, n_sims = n_sims, failed_file = failed_file,
       path_prefix = path_prefix)
}

#' Array indices of the grid rows matching a set of column values
#'
#' The safe way to run part of a study. Filtering the grid itself renumbers every
#' job (bug D); this selects row numbers and leaves the numbering alone.
#'
#'   idx <- grid_indices(study, method = "complete_data")
#'   cat(idx, file = "jobscripts/subset_ids.txt", sep = "\n")
#'
#' @param ... column = value(s) pairs to match
grid_indices <- function(study, ...) {
  conds <- list(...)
  keep <- rep(TRUE, nrow(study$grid))
  for (nm in names(conds)) {
    if (!nm %in% names(study$grid)) stop("no column '", nm, "' in the grid")
    keep <- keep & study$grid[[nm]] %in% conds[[nm]]
  }
  which(keep)
}

#' Parameter combinations, without the run index
combos <- function(study) {
  unique(study$grid[, study$path_cols, drop = FALSE])
}

#' Directory holding the per-run files for one parameter combination
combo_dir <- function(study, combo) {
  parts <- vapply(study$path_cols, function(cl) {
    prefix <- if (cl %in% names(study$path_prefix)) study$path_prefix[[cl]] else ""
    paste0(prefix, as.character(combo[[cl]]))
  }, character(1))
  do.call(file.path, c(list(study$res_path), as.list(parts)))
}

RES_PATTERN <- "^res_sim_\\d+\\.RDS$"

run_numbers <- function(files) {
  as.integer(gsub(".*res_sim_(\\d+)\\.RDS$", "\\1", files))
}

#' Which runs are missing, and the array indices that would produce them
#'
#' Writes the indices to study$failed_file so they can be resubmitted directly.
#' Indices are row numbers of study$grid, which is why the grid must never be
#' filtered after construction.
#'
#' Also points <prefix>_rerun.sh at the list it just wrote - see
#' update_rerun_script(). Keeping the two in step used to be a manual edit the
#' rerun scripts asked for in a comment, and nothing checked it: -J below the
#' line count silently drops the tail of the list, and -J above it feeds the
#' analysis script an empty index.
#'
#' @param write write the index list to study$failed_file
#' @param update_rerun also rewrite the rerun jobscript. Only consulted when
#'   `write` is TRUE, so check_failed(study, write = FALSE) touches nothing.
#' @param ... passed to update_rerun_script(), e.g. mem_factor, time_factor
check_failed <- function(study, write = TRUE, update_rerun = TRUE, ...) {

  cmb <- combos(study)

  failed <- bind_rows(lapply(seq_len(nrow(cmb)), function(i) {
    combo <- cmb[i, , drop = FALSE]
    files <- list.files(combo_dir(study, combo), pattern = RES_PATTERN,
                        full.names = TRUE)
    if (length(files) >= study$n_sims) return(NULL)
    missing_runs <- setdiff(seq_len(study$n_sims), run_numbers(files))
    if (!length(missing_runs)) return(NULL)
    combo[rep(1, length(missing_runs)), , drop = FALSE] %>%
      mutate(run = missing_runs)
  }))

  if (nrow(failed) == 0) {
    print(paste("All simulations complete! Go ahead and collect up the",
                "results (rerun jobscript left as it is - nothing to submit,",
                "and -J 1-0 is not a legal array)"))
    return(invisible(integer(0)))
  }

  key_cols <- c(study$path_cols, "run")
  grid_key <- do.call(paste, c(study$grid[, key_cols], sep = "\r"))
  failed_key <- do.call(paste, c(failed[, key_cols], sep = "\r"))
  failed_idx <- which(grid_key %in% failed_key)

  if (write && !is.null(study$failed_file)) {
    cat(failed_idx, file = study$failed_file, sep = "\n")
    print(paste0("failed runs found (", nrow(failed),
                 ") saved to jobscripts folder"))
    if (update_rerun) update_rerun_script(study, length(failed_idx), ...)
  }

  invisible(failed_idx)
}

##### jobscript editing #######################################################
# check_failed() writes failed_ids.txt and then rewrites the study's
# <prefix>_rerun.sh to match it. Doing that by hand is what the "Set -J to
# 1-<number of lines in failed_ids.txt> before submitting" comments in the
# rerun scripts were asking for, and it is easy to get wrong in a way nothing
# catches: -J too small silently drops the tail of failed_ids.txt, and -J too
# large makes `sed -n "${PBS_ARRAY_INDEX}p"` return an empty string, so the
# analysis script is called with an empty index argument.

#' Read a jobscript without disturbing its line endings
#'
#' Every jobscript in this repo is CRLF, four of the rerun scripts have no
#' trailing newline, and there is no .gitattributes. readLines()/writeLines()
#' - the idiom the *_profile_summary.R scripts use - strips the \r from only
#' the lines it edits when it runs on the HPC login node, and always appends a
#' final newline, so a three-line edit shows up as a whole-file diff. Round-trip
#' the bytes instead and put back exactly what was there.
read_script <- function(path) {
  txt <- rawToChar(readBin(path, "raw", file.size(path)))
  list(lines    = strsplit(gsub("\r\n", "\n", txt, fixed = TRUE), "\n",
                           fixed = TRUE)[[1]],
       eol      = if (grepl("\r\n", txt, fixed = TRUE)) "\r\n" else "\n",
       final_nl = grepl("\n$", txt))
}

#' Write back what read_script() read, with its original endings
write_script <- function(sc, path) {
  txt <- paste0(paste(sc$lines, collapse = sc$eol),
                if (sc$final_nl) sc$eol else "")
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(charToRaw(txt), con)
}

#' Seconds in an HH:MM:SS walltime, or NULL if the line isn't one
parse_walltime <- function(line) {
  m <- regmatches(line, regexec("walltime=(\\d+):(\\d{2}):(\\d{2})", line))[[1]]
  if (!length(m)) return(NULL)
  n <- as.integer(m[-1])
  n[1] * 3600L + n[2] * 60L + n[3]
}

format_walltime <- function(secs) {
  secs <- round(secs)   # %d rejects a fractional double; a walltime is whole
  sprintf("%02d:%02d:%02d", secs %/% 3600, (secs %% 3600) %/% 60, secs %% 60)
}

#' ncpus / ompthreads / mem from a PBS select= line, or NULL if any is missing
#'
#' Each field is matched on its own so the order within select= doesn't matter.
#' mem must be a whole number of gb - every jobscript in the repo is written
#' that way, and silently mis-parsing a unit would produce a request that is
#' wrong by a factor of 1000 in whichever direction hurts more.
parse_select <- function(line) {
  one <- function(pat) {
    m <- regmatches(line, regexec(pat, line))[[1]]
    if (length(m)) as.integer(m[2]) else NA_integer_
  }
  out <- list(ncpus      = one("ncpus=(\\d+)"),
              ompthreads = one("ompthreads=(\\d+)"),
              mem        = one("mem=(\\d+)gb"))
  if (anyNA(unlist(out))) NULL else out
}

#' Point the rerun jobscript at the list check_failed() just wrote
#'
#' -J becomes 1-<n_failed>, and the resource request is lifted above the base
#' <prefix>_1.sh: a rerun exists because the first attempt failed, and the usual
#' reasons are the walltime and the memory limit, so resubmitting on identical
#' resources reruns the same failure.
#'
#' The bump is always computed from the base script and then max()ed against
#' what the rerun script already asks for. That makes repeat calls idempotent -
#' the second call finds its own output and changes nothing - while leaving hand
#' escalations alone: bin_ci_rerun.sh was raised to mem=15gb against a 5gb base,
#' and recomputing from the base alone would quietly drop it back to 6gb.
#'
#' Scripts are found by study$prefix rather than by globbing for *_rerun.sh,
#' because confidence_intervals/optimal_sf/jobscripts holds two studies and a
#' glob there matches both.
#'
#' @param n_failed number of indices in study$failed_file
#' A directive line in the rerun script that is present but unreadable is left
#' exactly as it is, with a warning, because there is then no way to tell what
#' it currently asks for and so no way to promise not to lower it. The other
#' directive is still updated - the two are independent.
#'
#' @param cpu_add cores added to the base ncpus; rounded up to a whole core.
#'   ompthreads is deliberately not raised with it: the R-level parallelism is
#'   fixed by the trailing args on the Rscript line (workers x grf_threads =
#'   base ncpus), which this does not touch, so the extra core is headroom
#'   rather than a thread nobody spawns.
#' @param mem_factor,time_factor multipliers on the base mem and walltime. Both
#'   results are rounded up - PBS wants whole gb and whole seconds.
#' @param max_walltime ceiling on the *computed* walltime. Applied before the
#'   max() against the current value, so a deliberately longer setting already
#'   in the script is never cut back to it. A value that is not HH:MM:SS warns
#'   and disables the cap rather than stopping the run.
#' @param throttle the %N concurrency limit written into -J
update_rerun_script <- function(study, n_failed, cpu_add = 1, mem_factor = 1.2,
                                time_factor = 2, max_walltime = "72:00:00",
                                throttle = 100) {

  if (is.null(study$prefix)) {
    warning("study has no prefix - cannot find its rerun script")
    return(invisible(NULL))
  }

  dir <- dirname(study$failed_file)
  rerun_path <- file.path(dir, paste0(study$prefix, "_rerun.sh"))
  base_path  <- file.path(dir, paste0(study$prefix, "_1.sh"))

  if (!file.exists(rerun_path)) {
    warning("no rerun script at ", rerun_path, " - -J not updated")
    return(invisible(NULL))
  }

  sc <- read_script(rerun_path)
  sel_i <- grep("^#PBS -l select=", sc$lines)
  wt_i  <- grep("^#PBS -l walltime=", sc$lines)
  j_i   <- grep("^#PBS -J ", sc$lines)

  # -J first: it is the half that must be right for the job to run at all, and
  # it is still worth writing when the base script is missing or unparseable.
  if (!length(j_i)) {
    warning("no '#PBS -J' line in ", rerun_path)
  } else {
    sc$lines[j_i[1]] <- sprintf("#PBS -J 1-%d%%%d", n_failed, throttle)
  }

  bumped <- NULL
  if (!file.exists(base_path)) {
    warning("no base script at ", base_path, " - resources left as they are")
  } else {
    base_sc <- read_script(base_path)
    base_sel <- parse_select(grep("^#PBS -l select=", base_sc$lines,
                                  value = TRUE)[1])
    base_wt  <- parse_walltime(grep("^#PBS -l walltime=", base_sc$lines,
                                    value = TRUE)[1])
    cur_sel <- if (length(sel_i)) parse_select(sc$lines[sel_i[1]]) else NULL
    cur_wt  <- if (length(wt_i))  parse_walltime(sc$lines[wt_i[1]]) else NULL

    if (is.null(base_sel) || is.null(base_wt)) {
      warning("could not read the resource request in ", base_path,
              " - resources left as they are")
    } else {
      # A directive that is present but unreadable is the dangerous case: we
      # cannot tell what it currently asks for, so we cannot honour the promise
      # not to lower it. Leave that line alone and say so, rather than read a
      # missing floor as a floor of zero and quietly shrink the request -
      # mem=15000mb is the same ask as mem=15gb, and parse_select() rejecting
      # it must not cost 9gb. A line that is simply absent is different: there
      # is nothing to protect and nothing to write.
      sel_stuck <- length(sel_i) > 0 && is.null(cur_sel)
      wt_stuck  <- length(wt_i)  > 0 && is.null(cur_wt)
      if (sel_stuck) {
        warning("cannot read the current select= in ", basename(rerun_path),
                " - left as it is; write it as mem=<n>gb to let this update it")
      }
      if (wt_stuck) {
        warning("cannot read the current walltime= in ", basename(rerun_path),
                " - left as it is")
      }

      floor_sel <- if (is.null(cur_sel)) {
        list(ncpus = 0L, ompthreads = 0L, mem = 0L)
      } else cur_sel

      # ceiling(), not the raw product: PBS wants whole cores and whole seconds,
      # and sprintf's %d rejects a fractional double outright, so a non-integer
      # cpu_add or time_factor used to abort the call before anything was
      # written - taking the -J fix down with it.
      ncpus <- max(ceiling(base_sel$ncpus + cpu_add), floor_sel$ncpus)
      omp   <- max(base_sel$ompthreads,               floor_sel$ompthreads)
      mem   <- max(ceiling(base_sel$mem * mem_factor), floor_sel$mem)

      cap <- parse_walltime(paste0("walltime=", max_walltime))
      if (is.null(cap)) {
        warning("max_walltime '", max_walltime,
                "' is not HH:MM:SS - no cap applied")
        cap <- Inf   # comparison stays well-defined; the cap branch never fires
      }
      wt <- ceiling(base_wt * time_factor)
      if (wt > cap) {
        warning(study$name, ": walltime x", time_factor, " is ",
                format_walltime(wt), ", capped at ", max_walltime)
        wt <- cap
      }
      wt <- max(wt, if (is.null(cur_wt)) 0 else cur_wt)

      if (length(sel_i) && !sel_stuck) {
        sc$lines[sel_i[1]] <- sprintf(
          "#PBS -l select=1:ncpus=%d:ompthreads=%d:mem=%dgb", ncpus, omp, mem)
      }
      if (length(wt_i) && !wt_stuck) {
        sc$lines[wt_i[1]] <- sprintf("#PBS -l walltime=%s", format_walltime(wt))
      }
      # report what actually landed, not what was computed - a skipped line
      # would otherwise be announced as though it had been applied
      parts <- c(
        if (length(sel_i) && !sel_stuck)
          sprintf("ncpus=%d ompthreads=%d mem=%dgb", ncpus, omp, mem),
        if (length(wt_i) && !wt_stuck)
          sprintf("walltime=%s", format_walltime(wt)))
      bumped <- if (length(parts)) paste(parts, collapse = " ") else NULL
    }
  }

  write_script(sc, rerun_path)
  print(paste0(basename(rerun_path), ": -J 1-", n_failed, "%", throttle,
               if (is.null(bumped)) "" else paste0(", ", bumped)))

  invisible(rerun_path)
}

#' Read every per-run result file into a nested tibble
#'
#' One row per parameter combination, with a list-column of
#' list(run = , result = ) entries - the shape the *_metrics.R scripts unnest.
get_results <- function(study, workers = 2) {

  cmb <- combos(study)

  read_combo <- function(i) {
    combo <- cmb[i, , drop = FALSE]
    files <- list.files(combo_dir(study, combo), pattern = RES_PATTERN,
                        full.names = TRUE)
    if (length(files) == 0) return(NULL)

    temp <- map(files, function(f) {
      list(run = run_numbers(f), result = readRDS(f))
    })
    gc()
    bind_cols(as_tibble(combo), tibble(results = list(temp)))
  }

  if (workers > 1) {
    require(future); require(furrr)
    future::plan(future::multisession, workers = workers)
    on.exit(future::plan(future::sequential), add = TRUE)
    all_results <- furrr::future_map(seq_len(nrow(cmb)), read_combo)
  } else {
    all_results <- lapply(seq_len(nrow(cmb)), read_combo)
  }

  bind_rows(all_results[!vapply(all_results, is.null, logical(1))])
}
