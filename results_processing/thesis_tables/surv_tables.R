##########
# title: LaTeX tables for the thesis chapter - competing risks
##########
# The competing-risks study (competing_risk/), two kinds of table, each at the
# primary rho = 0.5 and the sensitivity rho = 0:
#   main       one arm per estimator family (the production arm), so the
#              families - and so the ways of handling the competing event -
#              compare across a row block
#   secondary  one table per family that ships in several arms (the
#              pseudo-value construction x fitting comparison, ADEMP.md "Aims"):
#              each arm's raw values, then the paired differences between arms
#              that differ by one factor
# Rows are scenario / censoring / arm, the columns CATE bias, RMSE and Pearson
# for Event 1 then Event 2, each cell "mean (MCSE)". The scenario, censoring
# and main-table arm labels (SURV_*) and the summary come from R/figures.R,
# shared with thesis_figures/surv.R; the secondary tables' arm labels are
# inline below. The layout is R/tables.R's surv_metrics_table().
#
# Writes to ../results/thesis_tables/, for <rho> in rho05, rho0:
#   surv_main_<rho>.tex
#   surv_<family>_<rho>.tex   family in pseudo_cf, pseudo_dr, rsf_dr, sl_t, sl_dr
# \input{} them into a document with
#   \usepackage{booktabs, longtable, pdflscape, array}
#
# Metric choices:
# - Bias is the CATE bias, `bias` (mean of est - true over units), as in
#   ss_tables.R. `ate_bias` is the same number unless a run has NA estimates.
# - RMSE, the per-run root mean squared error averaged over runs (so not the
#   square root of the mean MSE), in days.
# - Pearson. Scenario 1 is not null here (its CATE varies with X1 and X2), so
#   its correlation is a real measurement, from a surv_metrics.R run since
#   2026-10-06; an older metrics file has a placeholder 0 there and is refused.
#   Pearson is NA, and prints as a dash, where the arm's truth is constant:
#   the net RMST CATE of event 1 in scenarios 2 and 5 and of event 2 in 1 and 3
#   (ipw and csf_cs only).
# - "Combined" (RMSTc) is left out, as in surv_results.qmd.
#
# The arms are not all scored against the same truth (ADEMP.md, "Estimands"):
# ipw and csf_cs against the net (cause-specific) RMST CATE, csf_sh against the
# subdistribution RMST CATE (= -RMTL CATE), every pseudo-value arm against the
# RMTL CATE. The main table marks the net-RMST rows with a dagger.
#
# The paired differences: every arm of a run is fitted to the same data, so
# each metric is differenced per run and target, and its MCSE is
# sd(diff) / sqrt(pairs), as paired_rho_diff() does across rho. Runs missing
# from either arm drop out of the pairs.

library(here)
source(here("R", "figures.R"))
source(here("R", "tables.R"))

# paths
path <- here()
res_path <- file.path(dirname(path), "results", "competing_risk")
tab_path <- file.path(dirname(path), "results", "thesis_tables")
dir.create(tab_path, showWarnings = FALSE, recursive = TRUE)

# ---- labels ------------------------------------------------------------------
# scenario, censoring and production-arm labels are SURV_* in R/figures.R,
# shared with thesis_figures/surv.R

TARGETS <- c("Event 1", "Event 2")

RHOS <- list(
  list(out = "rho05", rho = 0.5, text = "$\\rho = 0.5$ (primary analysis)"),
  list(out = "rho0", rho = 0, text = "$\\rho = 0$ (sensitivity analysis)")
)

# the main table's rows: the production arm of each family, in row order, the
# net-RMST arms daggered
MAIN_ARMS <- SURV_ARM_LABELS
MAIN_ARMS[SURV_NET_ARMS] <- paste0(MAIN_ARMS[SURV_NET_ARMS], "$^\\dagger$")

# the secondary tables: each family's arms (raw rows, in order), then the
# contrasts (arm - ref), each between two arms that differ by one factor:
# fitting with the pseudo-values held at whole-sample, or the pseudo-values
# with the fitting held at single crossfit
FIT_LABEL <- "$\\Delta$ fitting: single CF $-$ OOB"
PV_LABEL <- "$\\Delta$ PV: crossfit $-$ whole"

FAMILIES <- list(
  list(
    out = "pseudo_cf",
    label = "Causal forest on pseudo-values",
    arms = c(pseudo_cf_whole_oob = "Whole PV, OOB",
             pseudo_cf_whole_scf = "Whole PV, single CF",
             pseudo_cf_cvps_scf = "Crossfit PV, single CF"),
    contrasts = tibble::tribble(
      ~arm,                  ~ref,                  ~label,
      "pseudo_cf_whole_scf", "pseudo_cf_whole_oob", FIT_LABEL,
      "pseudo_cf_cvps_scf",  "pseudo_cf_whole_scf", PV_LABEL
    ),
    note = ""
  ),
  list(
    out = "pseudo_dr",
    label = "DR-learner, random forests on pseudo-values",
    arms = c(pseudo_dr_whole_oob = "Whole PV, OOB",
             pseudo_dr_whole_scf = "Whole PV, single CF",
             pseudo_dr_cvps_scf = "Crossfit PV, single CF"),
    contrasts = tibble::tribble(
      ~arm,                  ~ref,                  ~label,
      "pseudo_dr_whole_scf", "pseudo_dr_whole_oob", FIT_LABEL,
      "pseudo_dr_cvps_scf",  "pseudo_dr_whole_scf", PV_LABEL
    ),
    note = paste0(
      " Crossfit pseudo-values only train the nuisance regressions: the DR ",
      "correction term uses the whole-sample pseudo-values in every arm."
    )
  ),
  list(
    out = "rsf_dr",
    label = "DR-learner, competing-risks random survival forest",
    arms = c(rsf_dr_oob = "OOB",
             rsf_dr_scf = "Single CF"),
    contrasts = tibble::tribble(
      ~arm,         ~ref,         ~label,
      "rsf_dr_scf", "rsf_dr_oob", "$\\Delta$ single CF $-$ OOB"
    ),
    note = paste0(
      " The outcome model is fitted to the observed times and causes, so ",
      "pseudo-values (whole-sample) enter the DR correction term only, and the ",
      "arms differ in fitting alone."
    )
  ),
  list(
    out = "sl_t",
    label = "T-learner, SuperLearner on pseudo-values",
    arms = c(sl_t_whole = "Whole PV",
             sl_t_cvps = "Crossfit PV",
             sl_t_split = "Split PV"),
    contrasts = tibble::tribble(
      ~arm,         ~ref,         ~label,
      "sl_t_cvps",  "sl_t_whole", PV_LABEL,
      "sl_t_split", "sl_t_cvps",  "$\\Delta$ split $-$ crossfit PV"
    ),
    note = paste0(
      " All arms are single crossfit. Split: the pseudo-values recomputed on, ",
      "and the model trained on, $V - 2$ folds rather than $V - 1$."
    )
  ),
  list(
    out = "sl_dr",
    label = "DR-learner, SuperLearner on pseudo-values",
    arms = c(sl_dr_whole = "Whole PV",
             sl_dr_cvps = "Crossfit PV"),
    contrasts = tibble::tribble(
      ~arm,         ~ref,          ~label,
      "sl_dr_cvps", "sl_dr_whole", PV_LABEL
    ),
    note = paste0(
      " Both arms are single crossfit. Crossfit pseudo-values only train the ",
      "nuisance regressions: the DR correction term uses the whole-sample ",
      "pseudo-values in both."
    )
  )
)

# the columns of each event, left to right. Drop a row here to drop a column.
# Pearson gets a third decimal in the secondary tables: the paired differences
# are hundredths.
table_cols <- function(digits_corr) {
  tibble::tribble(
    ~stem,  ~header,   ~digits,
    "bias", "Bias",    3,
    "rmse", "RMSE",    3,
    "corr", "Pearson", digits_corr
  )
}

CENSORING_NOTE <- paste0(
  "Censoring yes: uniform censoring on $(1, 180)$ as well as administrative ",
  "censoring at 180; no: administrative only."
)

# ---- data --------------------------------------------------------------------

metrics <- readRDS(file.path(res_path, "surv_metrics.RDS"))

if (!"rho" %in% names(metrics)) {
  stop("surv_metrics.RDS has no rho column: it predates the 2026-10-01 grid. ",
       "Rerun surv_collect.R and surv_metrics.R on the current results.",
       call. = FALSE)
}
s1_corr <- metrics$corr[metrics$scenario == 1]
if (length(s1_corr) > 0 && all(s1_corr %in% 0)) {
  stop("scenario 1's Pearson correlations are all 0, cate_metrics()'s ",
       "null-scenario placeholder: surv_metrics.RDS predates the 2026-10-06 ",
       "fix. Rerun surv_metrics.R.", call. = FALSE)
}

metrics <- metrics %>%
  filter(target %in% TARGETS)

stems <- table_cols(2)$stem
group_cols <- c("scenario", "censoring", "arm", "target")
run_cols <- c("scenario", "censoring", "framework", "target")

# scenario, censoring, target and arm to display factors; `arms` names the
# frameworks to keep and their row labels, in order
label_rows <- function(df, arms) {
  df %>%
    filter(framework %in% names(arms)) %>%
    mutate(
      scenario = label_factor(scenario, SURV_SCENARIO_LABELS),
      censoring = label_factor(censoring, SURV_CENSORING_LABELS),
      target = factor(target, levels = TARGETS),
      arm = factor(arms[framework], levels = unname(arms))
    )
}

#' Per-run arm - ref differences of `stems`, one row per run pair
paired_arm_diff <- function(metrics, arm, ref, stems,
                            keys = c("rho", "scenario", "n", "censoring",
                                     "run", "target")) {
  ref_runs <- metrics %>%
    filter(framework == ref) %>%
    select(all_of(keys), all_of(stems))
  out <- metrics %>%
    filter(framework == arm) %>%
    select(all_of(keys), all_of(stems)) %>%
    inner_join(ref_runs, by = keys, suffix = c("", "_ref"))
  for (s in stems) out[[s]] <- out[[s]] - out[[paste0(s, "_ref")]]
  select(out, all_of(keys), all_of(stems))
}

summarise_surv <- function(df) {
  summarise_metrics(df, group_cols, cols = setNames(stems, stems),
                    count_na = character())
}

write_table <- function(tex, out) {
  out_file <- file.path(tab_path, paste0(out, ".tex"))
  writeLines(tex, out_file)
  message("wrote ", out_file)
}

# ---- tables ------------------------------------------------------------------

for (r in RHOS) {
  m <- filter(metrics, rho == r$rho)
  if (nrow(m) == 0) {
    message("no rho = ", r$rho, " metrics - ", r$out, " tables skipped")
    next
  }

  # main: one arm per family
  main <- label_rows(m, MAIN_ARMS)
  out <- paste0("surv_main_", r$out)
  main %>%
    summarise_surv() %>%
    surv_metrics_table(
      table_cols(2),
      caption.short = paste0("CATE estimation by estimator family, ",
                             "competing-risks study, ", r$text),
      caption = paste0(
        "CATE estimation by estimator family, competing-risks study, $n = 500$, ",
        r$text, ": mean (Monte Carlo SE) over ",
        runs_per_cell(main, run_cols), " runs per cell. Event 1, event 2: ",
        "the CATE on that event's restricted mean time to the horizon of 28 ",
        "days. Rows marked $\\dagger$ are scored against the net ",
        "(cause-specific) RMST CATE, CSF subdistribution against the ",
        "subdistribution RMST CATE ($-$RMTL), and the others against the RMTL ",
        "CATE, so the sign of the bias and the scale of the errors differ ",
        "between those three groups. Bias is the CATE bias, the mean of the ",
        "estimated minus true CATE over units. ", CENSORING_NOTE,
        " --- : not defined (correlation with a constant true CATE)."
      ),
      label = out
    ) %>%
    write_table(out)

  # secondary: per family, the raw arms then the paired contrasts
  for (f in FAMILIES) {
    raw <- label_rows(m, f$arms)
    if (nrow(raw) == 0) {
      message(f$out, ": no ", r$out, " metrics - table skipped")
      next
    }

    diffs <- bind_rows(lapply(seq_len(nrow(f$contrasts)), function(i) {
      paired_arm_diff(m, f$contrasts$arm[i], f$contrasts$ref[i], stems) %>%
        mutate(framework = f$contrasts$label[i])
    }))
    diffs <- label_rows(diffs, setNames(f$contrasts$label, f$contrasts$label))

    row_labels <- c(unname(f$arms), f$contrasts$label)
    summary <- bind_rows(summarise_surv(raw), summarise_surv(diffs)) %>%
      mutate(arm = factor(as.character(arm), levels = row_labels))

    out <- paste0("surv_", f$out, "_", r$out)
    surv_metrics_table(
      summary,
      table_cols(3),
      caption.short = paste0(f$label, ", competing-risks study, ", r$text),
      caption = paste0(
        f$label, ", competing-risks study, $n = 500$, ", r$text, ". ",
        "The arm rows are means (Monte Carlo SE) over ",
        runs_per_cell(raw, run_cols), " runs per cell. The $\\Delta$ rows are ",
        "per-run differences between two arms, mean (Monte Carlo SE) over ",
        runs_per_cell(diffs, run_cols), " run pairs: every arm of a run is ",
        "fitted to the same data, so the SE is that of the paired difference.",
        if (FIT_LABEL %in% f$contrasts$label) {
          paste0(" Fitting: single crossfit minus whole-sample out-of-bag, ",
                 "both with whole-sample pseudo-values (PV).")
        },
        if (PV_LABEL %in% f$contrasts$label) {
          paste0(" PV: leave-one-fold-out minus whole-sample pseudo-values, ",
                 "both single crossfit.")
        },
        f$note,
        " Every arm is scored against the RMTL CATE. Bias is the CATE bias, ",
        "the mean of the estimated minus true CATE over units. ",
        CENSORING_NOTE
      ),
      label = out
    ) %>%
      write_table(out)
  }
}
