##########
# title: create every PBS log directory the jobscripts write to
##########
# Rscript make_log_dirs.R            # create any that are missing
# Rscript make_log_dirs.R --dry-run  # just list what would be created
#
# PBS aborts a job whose -o/-e directory doesn't exist, and *logs*/ is
# gitignored, so none of them survive a fresh clone. Run this on the HPC
# login node after cloning/pulling, before the first qsub.
#
# The directories are read from the #PBS -o / -e lines of every
# */jobscripts/*.sh rather than listed here, so a new jobscript is picked up
# without editing this file. A relative path is resolved against the
# jobscripts directory, because that is where qsub is run from (the scripts
# cd "${PBS_O_WORKDIR}/.." and read "${PBS_O_WORKDIR}/failed_ids.txt").

library(here)

dry_run <- "--dry-run" %in% commandArgs(trailingOnly = TRUE)

jobscripts <- list.files(here(), pattern = "\\.sh$", recursive = TRUE, full.names = TRUE)
jobscripts <- jobscripts[basename(dirname(jobscripts)) == "jobscripts"]

log_dirs <- unlist(lapply(jobscripts, function(f) {
  directives <- grep("^#PBS\\s+-[oe]\\s+", readLines(f, warn = FALSE), value = TRUE)
  paths <- trimws(sub("^#PBS\\s+-[oe]\\s+", "", directives))
  # -o/-e can name a file rather than a directory; only a trailing / means
  # "write into this directory", which is the form every jobscript here uses
  paths <- paths[grepl("/$", paths)]
  ifelse(grepl("^(/|~)", paths), paths, file.path(dirname(f), paths))
}))
log_dirs <- sort(unique(normalizePath(log_dirs, winslash = "/", mustWork = FALSE)))

missing_dirs <- log_dirs[!dir.exists(log_dirs)]
root <- paste0(normalizePath(here(), winslash = "/"), "/")
rel <- function(p) ifelse(startsWith(p, root), substring(p, nchar(root) + 1), p)

message(length(log_dirs), " log directories referenced by ", length(jobscripts),
        " jobscripts; ", length(missing_dirs), " missing")

for (d in missing_dirs) {
  if (dry_run) {
    message("  would create ", rel(d))
  } else if (dir.create(d, recursive = TRUE)) {
    message("  created ", rel(d))
  } else {
    warning("could not create ", d)
  }
}
