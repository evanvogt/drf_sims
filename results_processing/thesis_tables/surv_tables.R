##########
# title: LaTeX tables for the thesis chapter - competing risks
##########
# Appendix tables for the competing-risks study at rho = 0.5 and rho = 0: a main
# table with each estimator family's production arm (each DR-learner followed
# by its semi-oracle twin), and one secondary table
# per family with several arms (raw values, then paired differences between
# arms that differ by one factor). Cells are "mean (MCSE)" of CATE bias, RMSE
# and Pearson, EOI then CE.
#
# - Arms target different truths (ADEMP.md, "Estimands"): ipw and csf_cs the net
#   RMST CATE (daggered), csf_sh the subdistribution RMST, the pseudo-value arms
#   the RMTL.
# - RMSE is the per-run RMSE averaged over runs, in days.
# - Pearson is NA (a dash) where the arm's truth is constant. Scenario 1 is not
#   null here, so a metrics file from before 2026-10-06 (placeholder 0) is refused.
# - Every arm fits the same data, so differences are per run, MCSE
#   sd(diff) / sqrt(pairs).
# - "Combined" (RMSTc) is left out.
#
# Also the true CATEs the arms are scored against (surv_truth.tex): population
# mean (SD) of the RMTL and net RMST CATE per scenario, at both rho. Computed
# from the DGM, not read from the results, so it is written even without a
# metrics file; the grid truth takes a minute or two.
#
# Writes to ../results/thesis_tables/: surv_truth.tex, surv_main_<rho>.tex and
# surv_<family>_<rho>.tex, <rho> in rho05, rho0. \input{} with
# \usepackage{booktabs, longtable, array}; surv_truth.tex also uses
# the thesis's glossaries entries cate, rmtl and rmst.

library(here)
source(here("R", "figures.R"))
source(here("R", "tables.R"))
# DGM_SEED, DGM_N_POP(_SURV), dgm_params(); sources competing_risk/surv_dgm.R
source(here("results_processing", "dgm_population.R"))

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
# net-RMST arms daggered, and each DR-learner followed by its semi-oracle twin
# (true propensity 0.5 in the pseudo-outcome). Labels must be distinct: they
# become factor levels.
SEMI_ARMS <- c(
  pseudo_dr_whole_oob = "pseudo_dr_whole_oob_semi",
  rsf_dr_oob = "rsf_dr_oob_semi",
  sl_dr_whole = "sl_dr_whole_semi"
)
SEMI_LABELS <- c(
  pseudo_dr_whole_oob_semi = "\\quad DR-RF semi-oracle",
  rsf_dr_oob_semi = "\\quad DR-RSF semi-oracle",
  sl_dr_whole_semi = "\\quad DR-SL semi-oracle"
)

MAIN_ARMS <- SURV_ARM_LABELS
MAIN_ARMS[SURV_NET_ARMS] <- paste0(MAIN_ARMS[SURV_NET_ARMS], "$^\\dagger$")
MAIN_ARMS <- unlist(unname(lapply(names(MAIN_ARMS), function(a) {
  c(MAIN_ARMS[a], if (a %in% names(SEMI_ARMS)) SEMI_LABELS[SEMI_ARMS[[a]]])
})))

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

write_table <- function(tex, out) {
  out_file <- file.path(tab_path, paste0(out, ".tex"))
  writeLines(tex, out_file)
  message("wrote ", out_file)
}

# ---- true CATEs ----------------------------------------------------------------
# surv_truth_grid() per scenario, interpolated to a copula draw per rho, as
# surv_dgm_check.R does. The draw at each rho uses the same seed, so the two
# are paired, and at rho = 0.5 it is dgm_covariates()'s first DGM_N_POP_SURV
# rows, so the RMTL columns match dgm_tables.R's dgm_surv.tex.

TRUTH_COLS <- c(tau_RMTL1 = "RMTL1_0", tau_RMTL2 = "RMTL2_0",
                tau_RMST1_cs = "RMST1_cs_0", tau_RMST2_cs = "RMST2_cs_0")

truth_covariates <- function(rho) {
  p <- dgm_params("continuous", 1)
  p$rho <- rho
  set.seed(DGM_SEED)
  cv <- correlated_covariates(DGM_N_POP, p)
  lapply(cv[c("X1", "X2", "X3")], `[`, seq_len(DGM_N_POP_SURV))
}

#' "mean (SD)" with $-$ for minus, as fmt_mcse(); "0" where the CATE is
#' identically 0 (no effect on that cause's hazard)
fmt_truth <- function(v) {
  if (all(v == 0)) "0" else fmt_mcse(mean(v), sd(v), 2)
}

oplan <- future::plan(future::multisession,
                      workers = min(6, max(1, parallelly::availableCores() - 1)))
truth_grids <- lapply(survival_scenario_params$scenario, function(s) {
  surv_truth_grid(survival_scenario_params[survival_scenario_params$scenario == s, ])
})
future::plan(oplan)

truth_cells <- bind_rows(lapply(RHOS, function(r) {
  cv <- truth_covariates(r$rho)
  bind_rows(lapply(seq_along(truth_grids), function(s) {
    tibble(
      scenario = s,
      rho = r$rho,
      !!!lapply(setNames(names(TRUTH_COLS), names(TRUTH_COLS)), function(cl) {
        fmt_truth(surv_truth_at_draws(truth_grids[[s]], cl, cv$X1, cv$X2, cv$X3))
      })
    )
  }))
})) %>%
  arrange(scenario)  # stable: RHOS order within a scenario

# control-arm levels, the same in every scenario; at the primary rho
cv <- truth_covariates(RHOS[[1]]$rho)
ctrl <- sapply(TRUTH_COLS, function(cl) {
  sprintf("%.1f", mean(surv_truth_at_draws(truth_grids[[1]], cl, cv$X1, cv$X2, cv$X3)))
})

# scenario printed on its first row only, a line space between scenarios
truth_rows <- paste0(
  ifelse(duplicated(truth_cells$scenario), "", truth_cells$scenario), " & ",
  truth_cells$rho, " & ",
  do.call(paste, c(truth_cells[names(TRUTH_COLS)], sep = " & ")),
  " \\\\",
  ifelse(duplicated(truth_cells$scenario) &
           truth_cells$scenario < max(truth_cells$scenario), " \\addlinespace", "")
)

c(
  "\\begin{table}[h]",
  "\\small",
  paste0(
    "\\caption[True CATEs in the competing-risks simulation study]{True ",
    "\\glspl{cate} in the competing-risks simulation study: mean (SD) over ",
    format(DGM_N_POP_SURV, big.mark = ",", scientific = FALSE), " covariate ",
    "draws per $\\rho$, in days, at $t^* = 28$. $\\rho$: the covariates' ",
    "correlation, 0.5 in the primary analysis and 0 in the sensitivity ",
    "analysis. \\gls{rmtl}: positive when treatment brings the event forward. ",
    "Net \\gls{rmst}: negative when treatment brings the event forward. ",
    "Control-arm values ($\\rho = ", RHOS[[1]]$rho, "$): \\gls{rmtl} ",
    ctrl[["tau_RMTL1"]], " (EOI) and ", ctrl[["tau_RMTL2"]], " (CE); ",
    "net \\gls{rmst} ", ctrl[["tau_RMST1_cs"]], " and ", ctrl[["tau_RMST2_cs"]],
    ".}"
  ),
  "\\label{tab:surv_truth}",
  "\\begin{tabular}{@{}llllll@{}}",
  "\\toprule",
  paste0(" & & \\multicolumn{2}{l}{\\gls{rmtl}} & \\multicolumn{2}{l}{Net ",
         "\\gls{rmst}} \\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){5-6}"),
  "Scenario & $\\rho$ & EOI & CE & EOI & CE \\\\ \\midrule",
  truth_rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
) %>%
  paste(collapse = "\n") %>%
  write_table("surv_truth")

# ---- data --------------------------------------------------------------------

metrics <- readRDS(file.path(res_path, "surv_metrics.RDS"))

if (!"rho" %in% names(metrics)) {
  stop("surv_metrics.RDS has no rho column: it predates the 2026-10-01 grid. ",
       "Rerun study.R collect and surv_metrics.R on the current results.",
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
      target = factor(target, levels = TARGETS, labels = c("EOI", "CE")),
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
        runs_per_cell(main, run_cols), " runs per cell. EOI, CE: ",
        "the CATE on that event's restricted mean time to the horizon of 28 ",
        "days. Rows marked $\\dagger$ are scored against the net ",
        "(cause-specific) RMST CATE, CSF subdistribution against the ",
        "subdistribution RMST CATE ($-$RMTL), and the others against the RMTL ",
        "CATE, so the sign of the bias and the scale of the errors differ ",
        "between those three groups. Semi-oracle rows: the DR-learner above ",
        "with the true propensity 0.5 in place of the estimated one in the ",
        "pseudo-outcome (same outcome models), scored against the RMTL CATE. ",
        "Bias is the CATE bias, the mean of the ",
        "estimated minus true CATE over units. ", CENSORING_NOTE,
        " --- : not defined (correlation with a constant true CATE)."
      ),
      label = out,
      landscape = FALSE
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
      label = out,
      landscape = FALSE
    ) %>%
      write_table(out)
  }
}
