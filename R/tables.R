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

#' Booktabs longtable of summarised metrics, one row per scenario x n x model
#'
#' Rows are grouped under a scenario header (pack_rows), n is printed on the
#' first row of its block, and blocks are separated by \addlinespace. Every
#' metric column is "mean (MCSE)".
#'
#' @param summary output of summarise_metrics(), with `scenario`, `n` and
#'   `model` as factors in display order (apply_labels()), and
#'   mean_<stem>/mcse_<stem> for every stem in `cols`
#' @param cols data frame with one row per metric column: `stem`, `header`
#'   (LaTeX, not escaped), `group` (the spanning header above it: plain text,
#'   add_header_above() eats backslashes) and `digits`
#' @param caption,label passed to kbl(); `label` gets knitr's "tab:" prefix
#' @param caption.short the list-of-tables entry, so the long caption (which
#'   explains the cells) stays out of it
#' @param font_size in pt
#' @param tabcolsep inter-column padding, narrower than LaTeX's 6pt default so
#'   the wide tables fit a landscape page. Set inside the table's own group, so
#'   it does not leak into the rest of the document
#' @return the LaTeX source, as a single string
ss_latex_table <- function(summary, cols, caption, label, caption.short = "",
                           font_size = 8, tabcolsep = "3pt") {
  body <- summary %>% arrange(scenario, n, model)

  cells <- lapply(seq_len(nrow(cols)), function(i) {
    fmt_mcse(body[[paste0("mean_", cols$stem[i])]],
             body[[paste0("mcse_", cols$stem[i])]],
             cols$digits[i])
  })
  names(cells) <- cols$stem

  block <- paste(body$scenario, body$n)
  tab <- bind_cols(
    tibble(n = ifelse(duplicated(block), "", as.character(body$n)),
           model = as.character(body$model)),
    as_tibble(cells)
  )

  # a line space after the last row of each n block, except the scenario's
  # last (pack_rows puts its own space before the next scenario's header).
  # Added with row_spec() after pack_rows(), which drops kbl()'s `linesep`.
  new_block <- c(block[-1] != block[-length(block)], FALSE)
  new_scenario <- c(body$scenario[-1] != body$scenario[-length(block)], FALSE)
  space_after <- which(new_block & !new_scenario)

  groups <- rle(cols$group)
  header_above <- c(2, groups$lengths)
  names(header_above) <- c(" ", groups$values)

  scenario_index <- table(droplevels(body$scenario))

  out <- kbl(
    tab,
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    escape = FALSE,
    col.names = c("$n$", "Model", cols$header),
    align = c("l", "l", rep("r", nrow(cols))),
    linesep = "",
    caption = caption,
    caption.short = caption.short,
    label = label
  ) %>%
    add_header_above(header_above, escape = FALSE) %>%
    pack_rows(index = scenario_index, escape = FALSE) %>%
    row_spec(space_after, extra_latex_after = "\\addlinespace") %>%
    # later pages get "Table x: (continued)", not the whole caption again
    kable_styling(latex_options = "repeat_header", repeat_header_method = "replace",
                  font_size = font_size) %>%
    landscape()

  paste0("{\\setlength{\\tabcolsep}{", tabcolsep, "}\n", out, "\n}\n")
}
