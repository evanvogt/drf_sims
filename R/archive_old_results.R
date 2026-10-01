##########
# title: archive every study's pre-2026-09-26 results, one .tar per study
##########
# Everything under ../results predates bug O (the continuous baseline and bW
# calibration), the risk-difference binary DGM and the 2026-09-26 scenario
# renumbering (R/dgm_scenarios.R). None of it is what the current code would
# produce, and its scenario_<k>/ directories use the old numbers. Re-running
# into those trees would overwrite some files - missing/continuous's new
# scenario 2 lands on the old scenario 2's paths - and leave others behind to
# be collected under the wrong scenario, so every tree is archived first and
# every study re-runs from empty.
#
# Each study's res_path becomes one .tar under results/_archive/pre_2026-09-26/,
# named after its path with "/" as "__". Uncompressed, because every
# res_sim_*.RDS is gzip-compressed already (saveRDS()'s default): compressing
# again would save little and cost hours. What the .tar buys is one file per
# study instead of tens of thousands. A tree is deleted only once its .tar lists
# exactly the tree's files; if not, the .tar is removed and the tree kept.
#
#   Rscript R/archive_old_results.R                  # dry run - what would happen
#   qsub R/jobscripts/archive_old_results.sh         # the real run, on the cluster
#   Rscript R/archive_old_results.R --apply          # the same, interactively
#   Rscript R/archive_old_results.R --root <dir>     # a results/ other than
#                                                    # ../results (for testing)
#
# Restore one tree:  tar -xf _archive/pre_2026-09-26/<name>.tar  (from results/)
#
# Which trees: the res_path of every study in R/study_registry.R, plus
# model_evaluation's strategies and split trees, which share me_config.R. A
# res_path inside another (optimal_sf's sf_calibration/, inside the CI studies')
# goes into its parent's .tar. Left alone: competing_risk/, whose DGM is its own
# (competing_risk/surv_dgm.R) and was unchanged by the 2026-09-26 change (its
# own 2026-10-01 retune is archived separately, by hand - competing_risk/
# README.md), sample_size/correlated/'s two studies, which were first run after
# 2026-09-26 so have nothing old to archive (and archiving them would tar new
# results), and anything no study config names -
# the figure directories, and ../collected_metrics/, which is outside results/.
#
# Safe to re-run: a tree whose .tar already exists is skipped, and every
# archived tree adds a row to _archive/pre_2026-09-26/MANIFEST.csv.

suppressPackageStartupMessages(library(here))
source(here("R", "pipeline.R"))
source(here("R", "study_registry.R"))

ARCHIVE_DIR <- file.path("_archive", "pre_2026-09-26")
UNAFFECTED <- c("competing_risk", "correlated/continuous", "correlated/binary")

args <- commandArgs(trailingOnly = TRUE)
apply <- "--apply" %in% args
default_root <- file.path(dirname(here()), "results")
results_root <- if ("--root" %in% args) {
  args[match("--root", args) + 1]
} else {
  default_root
}
if (is.na(results_root) || !dir.exists(results_root)) {
  stop("no results directory at ", results_root, call. = FALSE)
}

tar_bin <- "tar"
if (!nzchar(Sys.which(tar_bin))) stop("tar not found on the PATH", call. = FALSE)

# ---- which trees ------------------------------------------------------------

studies <- study_registry[!study_registry$study_name %in% UNAFFECTED, ]
configs <- rbind(
  studies[, c("study_name", "config_path", "config_var")],
  data.frame(study_name = c("model_evaluation (strategies)",
                            "model_evaluation (split)"),
             config_path = "model_evaluation/me_config.R",
             config_var = c("study_strat", "study_split"))
)

# every config builds its res_path as file.path(dirname(here()), "results", ...);
# the part after that is the tree's path inside results_root
rel_paths <- vapply(seq_len(nrow(configs)), function(i) {
  res_path <- load_study(configs$config_path[i], configs$config_var[i])$res_path
  prefix <- paste0(default_root, "/")
  if (!startsWith(res_path, prefix)) {
    stop(configs$study_name[i], "'s res_path is not under ", default_root, ": ",
         res_path, call. = FALSE)
  }
  substring(res_path, nchar(prefix) + 1)
}, character(1))

trees <- data.frame(study = configs$study_name, rel = rel_paths)
trees <- trees[!duplicated(trees$rel), ]
nested <- vapply(trees$rel, function(r) {
  any(startsWith(r, paste0(setdiff(trees$rel, r), "/")))
}, logical(1))
trees <- trees[!nested, ]
trees$archive <- file.path(ARCHIVE_DIR, paste0(gsub("/", "__", trees$rel), ".tar"))

# ---- archive ------------------------------------------------------------------

# tar is run from results_root with relative paths throughout: GNU tar reads a
# "C:/..." archive path as host:path, and relative paths keep the member names
# the trees' paths inside results/, so an extract there restores them in place
owd <- setwd(results_root)
on.exit(setwd(owd), add = TRUE)

tree_files <- function(rel) {
  sort(file.path(rel, list.files(rel, recursive = TRUE, all.files = TRUE)))
}

#' Files in a .tar, without its directory entries
tar_files <- function(archive) {
  out <- system2(tar_bin, c("-tf", shQuote(archive)), stdout = TRUE)
  if (!is.null(attr(out, "status"))) return(NULL)
  sort(sub("^\\./", "", out[!endsWith(out, "/")]))
}

archive_tree <- function(rel, archive) {
  files <- tree_files(rel)
  status <- system2(tar_bin, c("-cf", shQuote(archive), shQuote(rel)))
  listed <- if (status == 0) tar_files(archive) else NULL
  if (!identical(listed, files)) {
    unlink(archive)
    stop("the .tar of ", rel, " does not list exactly its ", length(files),
         " files (tar exit status ", status, ") - removed it, left ", rel,
         " in place", call. = FALSE)
  }
  unlink(rel, recursive = TRUE)
  if (dir.exists(rel)) {
    stop("archived ", rel, " to ", archive, " but could not delete it",
         call. = FALSE)
  }
  length(files)
}

cat(if (apply) "ARCHIVING" else "DRY RUN (pass --apply to archive)",
    "under", normalizePath(results_root, winslash = "/"), "\n\n")

manifest <- file.path(ARCHIVE_DIR, "MANIFEST.csv")
for (i in seq_len(nrow(trees))) {
  rel <- trees$rel[i]
  archive <- trees$archive[i]
  label <- sprintf("%-38s %s", trees$study[i], rel)

  if (file.exists(archive)) {
    if (dir.exists(rel)) {
      stop(archive, " exists but so does ", rel, " - an earlier run stopped ",
           "part-way. Check the .tar against the tree by hand before going on.",
           call. = FALSE)
    }
    cat(label, ": already archived\n", sep = "")
    next
  }
  if (!dir.exists(rel)) {
    cat(label, ": nothing there\n", sep = "")
    next
  }

  files <- tree_files(rel)
  bytes <- sum(file.size(files))
  cat(sprintf("%s: %d files, %.2f GB -> %s\n", label, length(files),
              bytes / 1e9, archive))
  if (!apply) next

  dir.create(ARCHIVE_DIR, recursive = TRUE, showWarnings = FALSE)
  n <- archive_tree(rel, archive)
  row <- data.frame(study = trees$study[i], from = rel, archive = archive,
                    n_files = n, bytes = bytes,
                    archived_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  write.table(row, manifest, sep = ",", row.names = FALSE,
              col.names = !file.exists(manifest), append = file.exists(manifest))
  cat("  archived and deleted\n")
}

untouched <- setdiff(list.dirs(".", full.names = FALSE, recursive = FALSE),
                     c(sub("/.*", "", trees$rel), sub("/.*", "", ARCHIVE_DIR)))
if (length(untouched)) {
  cat("\nnot archived (no study config names them):",
      paste(sort(untouched), collapse = ", "), "\n")
}
