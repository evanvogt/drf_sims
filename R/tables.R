##########
# title: shared helpers for the thesis LaTeX tables
##########
# The table counterpart of R/figures.R: presentation only. Each table script
# keeps its own paths, filters and choice of columns, and summarises with
# R/figures.R's summarise_metrics(), so a table cell and a figure's error bar
# come from the same mean and MCSE.
#
# Output is booktabs LaTeX from kableExtra. The document needs
#   \usepackage{booktabs, longtable, pdflscape, array}

library(dplyr)
library(tidyr)
library(kableExtra)

#' Format a mean and its Monte Carlo SE into one cell, "0.012 (0.003)"
#'
#' Both are rounded to `digits` decimals. Minus signs are typeset as `$-$`
#' (a hyphen otherwise), and a value that rounds to zero loses its sign, so no
#' "-0.000". A missing mean (an undefined metric, e.g. correlation in the null
#' scenario) gives an em dash; a missing MCSE (one non-NA run) gives the mean
#' alone.
fmt_mcse <- function(mean, mcse, digits) {
  case_when(
    is.na(mean) ~ "---",
    is.na(mcse) ~ fmt_num(mean, digits),
    TRUE ~ paste0(fmt_num(mean, digits), " (", fmt_num(mcse, digits), ")")
  )
}

#' A mean alone, as fmt_mcse() prints it: `digits` decimals, `$-$` for minus,
#' no "-0.000", an em dash for a missing value. For tables that give the MCSEs
#' once in the caption (mcse_note()).
fmt_num <- function(x, digits) {
  r <- round(x, digits)
  r[!is.na(r) & r == 0] <- 0
  out <- sub("^-", "$-$", formatC(r, format = "f", digits = digits))
  ifelse(is.na(x), "---", out)
}

#' Runs per cell for a caption, read off the data rather than the design, so
#' failed runs show: "500", or "99--100" when cells differ
runs_per_cell <- function(metrics, group_cols = c("scenario", "n", "model")) {
  r <- range(count(metrics, across(all_of(group_cols)), name = "runs")$runs)
  if (r[1] == r[2]) as.character(r[1]) else paste0(r[1], "--", r[2])
}

#' Booktabs longtable with rows grouped two levels deep
#'
#' The outer level (`group`, e.g. scenario) is a bold header row (pack_rows).
#' The inner level (`block`, e.g. n) is `tab`'s first column, printed on the
#' first row of its block only, with \addlinespace between blocks.
#'
#' @param tab the displayed columns, sorted by group then block. The first two
#'   columns are row labels (left-aligned), the rest are cells (right-aligned)
#' @param col_names column headers, LaTeX, not escaped
#' @param group,block one value per row of `tab`; `group` a factor. `block`
#'   NULL: each row is its own block, so the first column is never blanked and
#'   no space separates the rows of a group
#' @param header_above optional spanning headers, as add_header_above() takes
#'   them, or a list of such vectors for several rows, the first nearest the
#'   column headers. Plain text only: add_header_above() strips one level of
#'   backslash even with escape = FALSE, so e.g. "\\%" comes out as a comment
#'   character
#' @param align column alignments, as kbl() takes them; NULL for the default,
#'   two left-aligned label columns then right-aligned cells
#' @param caption,label passed to kbl(); `label` gets knitr's "tab:" prefix
#' @param caption.short the list-of-tables entry, so the long caption (which
#'   explains the cells) stays out of it
#' @param font_size in pt
#' @param tabcolsep inter-column padding, narrower than LaTeX's 6pt default so
#'   wide tables fit. Set inside the table's own group, so it does not leak
#'   into the rest of the document. NULL leaves the document's own padding
#' @param landscape rotate onto a landscape page (pdflscape)
#' @return the LaTeX source, as a single string
grouped_longtable <- function(tab, col_names, group, block, header_above = NULL,
                              caption, label, caption.short = "", font_size = 8,
                              tabcolsep = "3pt", landscape = TRUE, align = NULL) {
  # a line space after the last row of each block, except the group's last
  # (pack_rows puts its own space before the next group's header). Added with
  # row_spec() after pack_rows(), which drops kbl()'s `linesep`.
  space_after <- integer()
  if (!is.null(block)) {
    key <- paste(group, block)
    tab[[1]] <- ifelse(duplicated(key), "", as.character(tab[[1]]))
    last <- length(key)
    new_block <- c(key[-1] != key[-last], FALSE)
    new_group <- c(group[-1] != group[-last], FALSE)
    space_after <- which(new_block & !new_group)
  }

  out <- kbl(
    tab,
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    escape = FALSE,
    col.names = col_names,
    align = if (is.null(align)) c("l", "l", rep("r", ncol(tab) - 2)) else align,
    linesep = "",
    caption = caption,
    caption.short = caption.short,
    label = label
  )
  if (!is.null(header_above)) {
    if (!is.list(header_above)) header_above <- list(header_above)
    for (h in header_above) out <- add_header_above(out, h, escape = FALSE)
  }
  # no indent under the group header: the bold header is enough, and the 1em
  # goes to the columns
  out <- pack_rows(out, index = table(droplevels(group)), escape = FALSE,
                   indent = FALSE)
  if (length(space_after) > 0) {
    out <- row_spec(out, space_after, extra_latex_after = "\\addlinespace")
  }
  out <- out %>%
    # later pages get "Table x: (continued)", not the whole caption again
    kable_styling(latex_options = "repeat_header", repeat_header_method = "replace",
                  font_size = font_size)
  if (landscape) out <- landscape(out)

  # kableExtra tags repeated rows with \vphantom{k} to tell them apart while it
  # edits them. Done editing, so drop the tags: in a cell they add a space.
  out <- gsub(" ?\\\\vphantom\\{[0-9]+\\}", "", out)
  # where row_spec() finds the row's text verbatim it leaves the row's own \\
  # after the inserted \addlinespace, an empty row that doubles the space
  out <- gsub("\\\\addlinespace\\\\\\\\", "\\\\addlinespace", out)

  if (is.null(tabcolsep)) return(paste0(out, "\n"))
  paste0("{\\setlength{\\tabcolsep}{", tabcolsep, "}\n", out, "\n}\n")
}

#' Metrics table: one row per scenario x n x model, one column per metric
#'
#' Rows are grouped under a scenario header, then by n. Every metric column is
#' "mean (MCSE)".
#'
#' @param summary output of summarise_metrics(), with `scenario`, `n` and
#'   `model` as factors in display order (apply_labels()), and
#'   mean_<stem>/mcse_<stem> for every stem in `cols`
#' @param cols data frame with one row per metric column: `stem`, `header`
#'   (LaTeX, not escaped) and `digits`; optionally `group`, the spanning header
#'   above it (plain text, see grouped_longtable())
#' @param ... passed to grouped_longtable(): caption, label, caption.short,
#'   font_size, tabcolsep, landscape
ss_latex_table <- function(summary, cols, ...) {
  body <- summary %>% arrange(scenario, n, model)

  cells <- lapply(seq_len(nrow(cols)), function(i) {
    fmt_mcse(body[[paste0("mean_", cols$stem[i])]],
             body[[paste0("mcse_", cols$stem[i])]],
             cols$digits[i])
  })
  names(cells) <- cols$stem

  tab <- bind_cols(
    tibble(n = as.character(body$n), model = as.character(body$model)),
    as_tibble(cells)
  )

  grouped_longtable(
    tab,
    col_names = c("$n$", "Model", cols$header),
    group = body$scenario,
    block = body$n,
    header_above = col_group_headers(cols),
    ...
  )
}

#' Spanning headers from a `cols` table's optional `group` column, as
#' add_header_above() takes them (two row-label columns first); NULL without one
col_group_headers <- function(cols) {
  if (!"group" %in% names(cols)) return(NULL)
  groups <- rle(cols$group)
  setNames(c(2, groups$lengths), c(" ", groups$values))
}

#' Wide table: rows grouped under a header, one column per (spanning header,
#' column) pair
#'
#' The appendix's compact layout. A cell is a single formatted string, so the
#' caller decides what goes in it (a mean alone, with the MCSEs given once in
#' the caption by mcse_note()). Rows are `group` headers (pack_rows) over
#' `row_cols`; columns are `col` (e.g. n) under up to two spanning levels.
#' A (row, column) pair with no cell prints as an em dash, and a row with no
#' defined cell at all is dropped.
#'
#' @param d tidy data frame: `group` (factor, display order), the `row_cols`,
#'   the `spans` and `col` (factors or characters, columns in level/sort
#'   order) and `cell` (character)
#' @param row_cols row label columns, left to right, optionally named by their
#'   header; their sort order is the row order
#' @param spans spanning header columns nearest the column headers first (so
#'   `c("metric", "mechanism")` puts mechanism on the top row). Plain text only
#'   (see grouped_longtable(): a backslash is eaten, so write two for one)
#' @param col the column header column
#' @param block TRUE: the first row label is a sub-block, printed once with
#'   space between blocks (grouped_longtable's `block`); FALSE: every row is its
#'   own block
#' @param ... passed to grouped_longtable(): caption, label, caption.short,
#'   font_size, tabcolsep, landscape
wide_table <- function(d, row_cols, spans = character(), col = "col",
                       block = FALSE, ...) {
  row_names <- if (is.null(names(row_cols))) row_cols else names(row_cols)
  row_cols <- unname(row_cols)
  key_cols <- c(rev(spans), col)
  specs <- d %>%
    distinct(across(all_of(key_cols))) %>%
    arrange(across(all_of(key_cols)))
  id <- function(x) do.call(paste, c(lapply(x[key_cols], as.character), sep = "\r"))
  specs$.id <- id(specs)
  d$.id <- id(d)

  wide <- d %>%
    select(group, all_of(row_cols), .id, cell) %>%
    pivot_wider(names_from = .id, values_from = cell) %>%
    arrange(group, across(all_of(row_cols)))
  cell_cols <- specs$.id
  wide[cell_cols] <- lapply(wide[cell_cols], function(x) coalesce(x, "---"))
  keep <- rowSums(sapply(wide[cell_cols], function(x) x != "---")) > 0
  wide <- wide[keep, ]

  tab <- bind_cols(
    wide %>% transmute(across(all_of(row_cols), as.character)),
    as_tibble(wide[cell_cols], .name_repair = "minimal")
  )
  header_above <- lapply(spans, function(s) {
    r <- rle(as.character(specs[[s]]))
    setNames(c(length(row_cols), r$lengths), c(" ", r$values))
  })

  grouped_longtable(
    tab,
    col_names = c(row_names, as.character(specs[[col]])),
    group = wide$group,
    block = if (block) as.character(wide[[row_cols[1]]]) else NULL,
    header_above = if (length(header_above) > 0) header_above,
    align = c(rep("l", length(row_cols)), rep("r", length(cell_cols))),
    ...
  )
}

#' Largest Monte Carlo SE of each metric, as caption text
#'
#' Replaces an SE in every cell (Morris, White & Crowther 2019 accept MCSEs
#' reported in summary when small against the differences discussed).
#'
#' @param summary output of summarise_metrics() with mcse_<stem> columns
#' @param stems named character vector: names are the stems, values the
#'   metric's name in the sentence
#' @param digits decimals per stem, recycled
#' @return e.g. "Monte Carlo SEs were at most 0.007 (bias), 0.011 (MSE)."
mcse_note <- function(summary, stems, digits = 3) {
  digits <- rep_len(digits, length(stems))
  worst <- vapply(seq_along(stems), function(i) {
    x <- summary[[paste0("mcse_", names(stems)[i])]]
    # rounded up, so "at most" holds in the displayed digits
    formatC(ceiling(max(x, na.rm = TRUE) * 10^digits[i] - 1e-9) / 10^digits[i],
            format = "f", digits = digits[i])
  }, character(1))
  paste0("Monte Carlo SEs were at most ",
         paste0(worst, " (", unname(stems), ")", collapse = ", "), ".")
}

#' Largest binomial Monte Carlo SE of a rejection rate or 0/1 coverage, as
#' caption text
#'
#' sqrt(p(1 - p) / R) follows from the estimate, so it is at most 0.5 /
#' sqrt(R). `runs` is the number of runs (a vector for several cell sizes).
binom_mcse_note <- function(runs, what = "Monte Carlo SE") {
  runs <- sort(unique(runs))
  paste0(
    what, " $= \\sqrt{p(1-p)/R}$ for an estimate $p$ over $R$ runs: ",
    "at most ",
    paste0(formatC(0.5 / sqrt(runs), format = "f", digits = 3), " for $R = ",
           runs, "$", collapse = " and "),
    "."
  )
}

#' Competing-risks metrics table: rows scenario x censoring x arm, columns
#' event x metric
#'
#' ss_latex_table()'s layout with the competing-risks study's rows: grouped
#' under a scenario header, then by censoring. Each event (`target`) gets a
#' spanning header over its metric columns, so both events of one arm read
#' across a row. A cell with no row in `summary` prints as an em dash.
#'
#' @param summary output of summarise_metrics(), with `scenario`, `censoring`,
#'   `arm` and `target` as factors in display order, and mean_<stem>/mcse_<stem>
#'   for every stem in `cols`
#' @param cols as ss_latex_table(), without `group`: the spanning headers are
#'   the targets
#' @param ... passed to grouped_longtable()
surv_metrics_table <- function(summary, cols, ...) {
  summary <- droplevels(summary)
  targets <- levels(summary$target)

  long <- bind_rows(lapply(seq_len(nrow(cols)), function(i) {
    summary %>%
      transmute(
        scenario, censoring, arm, target,
        stem = cols$stem[i],
        cell = fmt_mcse(.data[[paste0("mean_", cols$stem[i])]],
                        .data[[paste0("mcse_", cols$stem[i])]],
                        cols$digits[i])
      )
  }))

  # event-major column order: every metric of the first event, then the next
  cell_cols <- paste(rep(targets, each = nrow(cols)), cols$stem, sep = "__")
  body <- long %>%
    mutate(key = paste(target, stem, sep = "__")) %>%
    select(-target, -stem) %>%
    pivot_wider(names_from = key, values_from = cell)
  for (cc in setdiff(cell_cols, names(body))) body[[cc]] <- NA_character_
  body <- body %>%
    mutate(across(all_of(cell_cols), ~ coalesce(.x, "---"))) %>%
    arrange(scenario, censoring, arm)

  tab <- body %>%
    transmute(censoring = as.character(censoring), arm = as.character(arm),
              across(all_of(cell_cols)))

  grouped_longtable(
    tab,
    col_names = c("Censoring", "Estimator", rep(cols$header, length(targets))),
    group = body$scenario,
    block = body$censoring,
    header_above = c(" " = 2,
                     setNames(rep(nrow(cols), length(targets)), targets)),
    ...
  )
}

#' Model-evaluation proxy table: rows scenario x n x score, columns fold regime
#' x machine
#'
#' One measure per table, so the 24 proxies of a (scenario, n) fold into four
#' rows: the score type down, the scoring nuisance's fold regime (a spanning
#' header) and machine across. Under them, one row per reference (random pick,
#' best fixed choice), which does not depend on the fold regime or machine and
#' so is printed once, in a cell spanning those columns.
#'
#' @param summary output of summarise_metrics() grouped by scenario, n,
#'   metric (the score type), fold and machine, all factors in display order,
#'   with mean_<stem>/mcse_<stem>
#' @param refs the references summarised the same way, grouped by scenario,
#'   n and ref (a factor, the row order), with `label`, the row's text
#' @param stem the measure
#' @param digits decimals for every cell
#' @param ... passed to grouped_longtable()
me_proxy_table <- function(summary, refs, stem, digits, ...) {
  mean_col <- paste0("mean_", stem)
  mcse_col <- paste0("mcse_", stem)
  summary <- droplevels(summary)
  folds <- levels(summary$fold)
  machines <- levels(summary$machine)
  cell_cols <- paste(rep(folds, each = length(machines)), machines, sep = "__")

  proxies <- summary %>%
    transmute(
      scenario, n,
      row = as.character(metric),
      order = as.integer(metric),
      key = paste(fold, machine, sep = "__"),
      cell = fmt_mcse(.data[[mean_col]], .data[[mcse_col]], digits)
    ) %>%
    pivot_wider(names_from = key, values_from = cell)
  for (cc in setdiff(cell_cols, names(proxies))) proxies[[cc]] <- NA_character_
  proxies <- proxies %>%
    mutate(across(all_of(cell_cols), ~ coalesce(.x, "---")))

  # kbl() needs every row to have all its cells, so a reference row's other
  # cells are placeholders, removed (with their "&") from the LaTeX afterwards
  placeholder <- "SPANNEDCELL"
  references <- refs %>%
    transmute(
      scenario, n,
      row = label,
      order = nlevels(summary$metric) + as.integer(ref),
      cell = fmt_mcse(.data[[mean_col]], .data[[mcse_col]], digits)
    )
  references[[cell_cols[1]]] <- paste0("\\multicolumn{", length(cell_cols),
                                       "}{c}{", references$cell, "}")
  for (cc in cell_cols[-1]) references[[cc]] <- placeholder
  references$cell <- NULL

  body <- bind_rows(proxies, references) %>%
    arrange(scenario, n, order)
  tab <- body %>%
    transmute(n = as.character(n), row, across(all_of(cell_cols)))

  out <- grouped_longtable(
    tab,
    col_names = c("$n$", "Score", rep(machines, length(folds))),
    group = body$scenario,
    block = body$n,
    header_above = c(" " = 2,
                     setNames(rep(length(machines), length(folds)), folds)),
    ...
  )
  gsub(paste0(" *& *", placeholder), "", out)
}

#' Model-evaluation design-axis table: rows scenario x n x axis x level, one
#' column per measure
#'
#' Each (scenario, n) is a bold header row, then each design axis (score,
#' propensity score, data splitting, surrogate learner) a block of its levels,
#' headed "Factor" as in the thesis.
#'
#' @param summary output of summarise_metrics() grouped by scenario, n, axis
#'   and level, all factors in display order, with mean_<stem>/mcse_<stem> for
#'   every stem in `cols`
#' @param cols as ss_latex_table()
#' @param ... passed to grouped_longtable()
me_axis_table <- function(summary, cols, ...) {
  body <- summary %>% arrange(scenario, n, axis, level)

  cells <- lapply(seq_len(nrow(cols)), function(i) {
    fmt_mcse(body[[paste0("mean_", cols$stem[i])]],
             body[[paste0("mcse_", cols$stem[i])]],
             cols$digits[i])
  })
  names(cells) <- cols$stem

  group <- paste0(body$scenario, ", $n$ = ", body$n)
  group <- factor(group, levels = unique(group))

  tab <- bind_cols(
    tibble(axis = as.character(body$axis), level = as.character(body$level)),
    as_tibble(cells)
  )

  grouped_longtable(
    tab,
    col_names = c("Factor", "Level", cols$header),
    group = group,
    block = body$axis,
    header_above = col_group_headers(cols),
    ...
  )
}
