##########
# title: LaTeX table for the thesis appendix - crossfitting comparison
##########
# The numbers behind thesis_figures/crossfitting.R's RMSE-ratio figure: per
# arm, the CATE RMSE on the training sample (the estimand) and on the
# independent test sample, by scenario, then the mean runtime per replicate
# and its multiple of the family's double-crossfit arm. Rows grouped by
# estimator family, double crossfitting first, the arm the main analyses use
# starred (CF_ARMS, R/figures.R).
#
# - RMSE is the per-run RMSE averaged over runs, cells "mean (MCSE)".
# - On the test sample a crossfit arm averages its V fold models' predictions,
#   which flatters it against a whole-sample arm (crossfitting/README.md).
# - Runtime is averaged over every run and scenario: it barely moves with the
#   scenario.
#
# Writes ../results/thesis_tables/cf_appendix.tex. \input{} with
# \usepackage{booktabs, longtable, pdflscape, array}.
#
#   cd results_processing/thesis_tables; Rscript cf_tables.R
#   Rscript cf_tables.R [results root] [table dir]   # mock smoke run

library(here)
source(here("R", "figures.R"))
source(here("R", "tables.R"))

# paths
args <- commandArgs(trailingOnly = TRUE)
res_root <- if (length(args) >= 1) args[1] else file.path(dirname(here()), "results")
tab_path <- if (length(args) >= 2) {
  args[2]
} else {
  file.path(dirname(here()), "results", "thesis_tables")
}
dir.create(tab_path, showWarnings = FALSE, recursive = TRUE)

SETS <- c(train = "Training sample ($n = 500$)", test = "Test sample ($n = 2000$)")
DIGITS_RMSE <- 3

# ---- data --------------------------------------------------------------------

metrics <- readRDS(file.path(res_root, "crossfitting", "cf_metrics.RDS")) %>%
  inner_join(CF_ARMS, by = c("arm", "family")) %>%
  mutate(
    scenario = label_factor(scenario, CF_SCENARIO_LABELS),
    family = label_factor(family, CF_FAMILY_LABELS),
    arm = factor(arm, levels = CF_ARMS$arm),
    set = factor(set, levels = names(SETS))
  )

rmse <- metrics %>%
  summarise_metrics(c("family", "arm", "set", "scenario"), cols = c(rmse = "rmse")) %>%
  mutate(cell = fmt_mcse(mean_rmse, mcse_rmse, DIGITS_RMSE)) %>%
  arrange(set, scenario) %>%
  select(family, arm, set, scenario, cell) %>%
  pivot_wider(names_from = c(set, scenario), values_from = cell)

# one row per arm-run in each set, so the training rows alone carry every run
runtime <- metrics %>%
  filter(set == "train") %>%
  group_by(family, arm, reference) %>%
  summarise(time = mean(time_total, na.rm = TRUE), .groups = "drop") %>%
  group_by(family) %>%
  mutate(rel_time = time / time[reference]) %>%
  ungroup() %>%
  transmute(family, arm,
            time = formatC(time, format = "f", digits = 1),
            rel_time = formatC(rel_time, format = "f", digits = 2))

body <- rmse %>%
  inner_join(runtime, by = c("family", "arm")) %>%
  left_join(CF_ARMS %>% mutate(arm = factor(arm, levels = CF_ARMS$arm)) %>%
              select(arm, nuisances, final, production),
            by = "arm") %>%
  arrange(family, arm)

scenarios <- levels(droplevels(metrics$scenario))
# "Single effects + interaction" over two lines, or its column is twice a cell
scenario_headers <- sub(" + ", "\\\\ + ", scenarios, fixed = TRUE)
scenario_headers <- ifelse(scenario_headers == scenarios, scenarios,
                           paste0("\\shortstack{", scenario_headers, "}"))
cell_cols <- setdiff(names(rmse), c("family", "arm"))

tab <- body %>%
  transmute(
    nuisances,
    final = ifelse(production, paste0(final, "$^\\ast$"), final),
    across(all_of(cell_cols)),
    time,
    rel_time
  )

# ---- table -------------------------------------------------------------------

out <- "cf_appendix"
tex <- grouped_longtable(
  tab,
  col_names = c("Nuisances", "Final model", rep(scenario_headers, length(SETS)),
                "Seconds", "$\\times$ double CF"),
  group = droplevels(body$family),
  block = NULL,
  header_above = c(" " = 2, setNames(rep(length(scenarios), length(SETS)), SETS),
                   "Runtime" = 2),
  caption.short = "CATE RMSE and runtime by crossfitting procedure",
  caption = paste0(
    "CATE RMSE and runtime by crossfitting procedure, continuous outcome, ",
    "$V = 10$ folds: mean (Monte Carlo SE) over ",
    runs_per_cell(metrics, c("scenario", "arm", "set")), " runs per cell. ",
    "Every arm of a run shares one fold assignment. Training sample: the ",
    "CATE of the $n = 500$ units the arm was fitted to, the estimand. Test ",
    "sample: an independent draw of 2000 from the same DGM, where a crossfit ",
    "arm averages the predictions of its $V$ fold models. Double CF: ",
    "nuisances fitted over all 45 fold pairs; single CF: leave one fold out; ",
    "OOB: one model on the whole sample, each unit predicted by the trees ",
    "that did not see it. Runtime: mean seconds per replicate on one cluster ",
    "core, over all scenarios, and as a multiple of the family's double-CF ",
    "arm. $^\\ast$ the procedure used in the main analyses."
  ),
  label = out
)

out_file <- file.path(tab_path, paste0(out, ".tex"))
writeLines(tex, out_file)
message("wrote ", out_file)
