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
  num <- function(x) {
    x <- round(x, digits)
    x[!is.na(x) & x == 0] <- 0
    sub("^-", "$-$", formatC(x, format = "f", digits = digits))
  }
  case_when(
    is.na(mean) ~ "---",
    is.na(mcse) ~ num(mean),
    TRUE ~ paste0(num(mean), " (", num(mcse), ")")
  )
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
#' @param group,block one value per row of `tab`; `group` a factor
#' @param header_above optional spanning headers, as add_header_above() takes
#'   them. Plain text only: add_header_above() strips one level of backslash
#'   even with escape = FALSE, so e.g. "\\%" comes out as a comment character
#' @param caption,label passed to kbl(); `label` gets knitr's "tab:" prefix
#' @param caption.short the list-of-tables entry, so the long caption (which
#'   explains the cells) stays out of it
#' @param font_size in pt
#' @param tabcolsep inter-column padding, narrower than LaTeX's 6pt default so
#'   wide tables fit. Set inside the table's own group, so it does not leak
#'   into the rest of the document
#' @param landscape rotate onto a landscape page (pdflscape)
#' @return the LaTeX source, as a single string
grouped_longtable <- function(tab, col_names, group, block, header_above = NULL,
                              caption, label, caption.short = "", font_size = 8,
                              tabcolsep = "3pt", landscape = TRUE) {
  key <- paste(group, block)
  tab[[1]] <- ifelse(duplicated(key), "", as.character(tab[[1]]))

  # a line space after the last row of each block, except the group's last
  # (pack_rows puts its own space before the next group's header). Added with
  # row_spec() after pack_rows(), which drops kbl()'s `linesep`.
  last <- length(key)
  new_block <- c(key[-1] != key[-last], FALSE)
  new_group <- c(group[-1] != group[-last], FALSE)
  space_after <- which(new_block & !new_group)

  out <- kbl(
    tab,
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    escape = FALSE,
    col.names = col_names,
    align = c("l", "l", rep("r", ncol(tab) - 2)),
    linesep = "",
    caption = caption,
    caption.short = caption.short,
    label = label
  )
  if (!is.null(header_above)) {
    out <- add_header_above(out, header_above, escape = FALSE)
  }
  out <- out %>%
    # no indent under the group header: the bold header is enough, and the
    # 1em goes to the columns
    pack_rows(index = table(droplevels(group)), escape = FALSE, indent = FALSE) %>%
    row_spec(space_after, extra_latex_after = "\\addlinespace") %>%
    # later pages get "Table x: (continued)", not the whole caption again
    kable_styling(latex_options = "repeat_header", repeat_header_method = "replace",
                  font_size = font_size)
  if (landscape) out <- landscape(out)

  # kableExtra tags repeated rows with \vphantom{k} to tell them apart while it
  # edits them. Done editing, so drop the tags: in a cell they add a space.
  out <- gsub(" ?\\\\vphantom\\{[0-9]+\\}", "", out)

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

  header_above <- NULL
  if ("group" %in% names(cols)) {
    groups <- rle(cols$group)
    header_above <- c(2, groups$lengths)
    names(header_above) <- c(" ", groups$values)
  }

  grouped_longtable(
    tab,
    col_names = c("$n$", "Model", cols$header),
    group = body$scenario,
    block = body$n,
    header_above = header_above,
    ...
  )
}

#' Test table: rows scenario x test x model, one column per sample size
#'
#' The same cells as a metrics table, turned so that a test's power reads
#' across n. A row whose cells are all undefined (a test that is never run for
#' that model, or never defined in that scenario) is dropped rather than
#' printed as a row of dashes.
#'
#' @param summary output of summarise_metrics(), with `scenario`, `n` and
#'   `model` as factors in display order, and mean_<stem>/mcse_<stem> for every
#'   stem in `tests`
#' @param tests named character vector: names are the stems, values the
#'   display labels, in row order
#' @param digits decimals for every cell
#' @param ... passed to grouped_longtable()
ss_test_table <- function(summary, tests, digits = 3, ...) {
  long <- bind_rows(lapply(names(tests), function(stem) {
    summary %>%
      transmute(
        scenario, n, model,
        test = tests[[stem]],
        cell = fmt_mcse(.data[[paste0("mean_", stem)]],
                        .data[[paste0("mcse_", stem)]], digits)
      )
  })) %>%
    mutate(test = factor(test, levels = unname(tests)))

  body <- long %>%
    arrange(n) %>%
    pivot_wider(names_from = n, values_from = cell, names_sort = TRUE) %>%
    filter(if_any(-c(scenario, model, test), ~ !is.na(.x) & .x != "---")) %>%
    mutate(across(-c(scenario, model, test), ~ coalesce(.x, "---"))) %>%
    arrange(scenario, test, model)

  n_cols <- setdiff(names(body), c("scenario", "model", "test"))
  tab <- body %>%
    transmute(test = as.character(test), model = as.character(model),
              across(all_of(n_cols)))

  grouped_longtable(
    tab,
    col_names = c("Test", "Model", n_cols),
    group = body$scenario,
    block = body$test,
    header_above = c(" " = 2, "Sample size" = length(n_cols)),
    ...
  )
}
