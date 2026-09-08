##########
# title: one-PNG-per-candidate comparison of poster panel candidates
##########
# One file per candidate encoding, under a subfolder per panel/variant, so
# the whole set can be flicked through image-by-image (e.g. arrow keys in a
# photo viewer) instead of squinted at as one crowded contact-sheet grid.
# See explore_options.qmd for the same candidates as browser tabs, and
# ../ICTMC_poster_planning.md for what each candidate is trying to settle.
# Candidate-building logic lives in R/candidates.R so neither this script nor
# the .qmd defines a plot twice.

library(here)
source(here("results_processing", "poster_figs", "R", "candidates.R"))

gallery_dir <- file.path(dirname(here()), "results", "ICTMC_figs", "gallery")

#' Save each candidate plot as its own PNG under gallery_dir/<subdir>/<name>.png
#'
#' @param plots named list of ggplot objects, as returned by
#'   candidate_plots_sample_size()/candidate_plots_missing()
#' @param subdir path (relative to gallery_dir) identifying the panel and
#'   variant, e.g. "sample_size/both" or "missing_data/by_model"
save_candidates <- function(plots, subdir) {
  out_dir <- file.path(gallery_dir, subdir)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  for (nm in names(plots)) {
    save_fig(paste0(nm, ".png"), path = out_dir, plot = plots[[nm]])
  }
}

# sample size - one subfolder per outcome-scope variant
save_candidates(candidate_plots_sample_size(outcome = "both"), "sample_size/both")
save_candidates(candidate_plots_sample_size(outcome = "binary"), "sample_size/binary")
save_candidates(candidate_plots_sample_size(outcome = "continuous"), "sample_size/continuous")

# missing data - averaged over model vs split by model
save_candidates(candidate_plots_missing(by_model = FALSE), "missing_data/averaged")
save_candidates(candidate_plots_missing(by_model = TRUE), "missing_data/by_model")
