# Regenerates README.html for every README.md that already has a rendered
# .html sibling, via the local `quarto` CLI. Local-machine doc maintenance
# only (not a cluster script) — keeps the two from drifting apart, which is
# what caused stale rendered READMEs during the 2026-08 README audit.
#
# Usage:
#   Rscript render_readmes.R                     # every README.md with an existing README.html
#   Rscript render_readmes.R path/to/README.md    # just that one (repeatable)

args <- commandArgs(trailingOnly = TRUE)

targets <- if (length(args) >= 1) {
  args
} else {
  all_md <- list.files(".", pattern = "^README\\.md$", recursive = TRUE, full.names = TRUE)
  all_md[file.exists(sub("\\.md$", ".html", all_md))]
}

if (length(targets) == 0) {
  stop("No README.md with an existing README.html sibling found.")
}

failures <- character(0)

for (f in targets) {
  message("Rendering ", f)
  status <- system2("quarto", c("render", shQuote(f)))
  if (status != 0) {
    warning("quarto render failed for ", f)
    failures <- c(failures, f)
  }
}

if (length(failures) > 0) {
  stop("quarto render failed for: ", paste(failures, collapse = ", "))
}

message("Rendered ", length(targets), " README(s) successfully.")
