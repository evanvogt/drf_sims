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
study_config <- function(name, res_path, grid, path_cols, n_sims = 100,
                         failed_file = NULL,
                         path_prefix = c(scenario = "scenario_")) {
  stopifnot("run" %in% names(grid))
  stopifnot(all(path_cols %in% names(grid)))
  list(name = name, res_path = res_path, grid = grid, path_cols = path_cols,
       n_sims = n_sims, failed_file = failed_file, path_prefix = path_prefix)
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
check_failed <- function(study, write = TRUE) {

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
    print("All simulations complete! Go ahead and collect up the results")
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
  }

  invisible(failed_idx)
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
