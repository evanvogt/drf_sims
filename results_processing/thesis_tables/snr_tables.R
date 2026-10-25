##########
# title: LaTeX table for the thesis chapter - binary SNR check, population truth
##########
# Appendix table for sample_size/snr/ (bin_snr_metrics.R): per scenario and HTE
# multiplier, the population truth behind the signal-to-noise axis. Reads only
# bin_snr_population.RDS; no simulation output.
#
# The competing-risk check's design table (hazard ratios, SNR) is typed in from
# competing_risk/snr/README.md: that study saves no population file.
#
# Writes to ../results/thesis_tables/:
#   bin_snr_population.tex
#   surv_snr_population.tex
# \input{} with \usepackage{booktabs, longtable, array}.
#
#   cd results_processing/thesis_tables; Rscript snr_tables.R
#   Rscript snr_tables.R [population dir] [table dir]   # mock smoke run

library(here)
source(here("R", "figures.R"))
source(here("R", "tables.R"))

# paths
path <- here()
args <- commandArgs(trailingOnly = TRUE)
res_path <- if (length(args) >= 1) args[1] else
  file.path(dirname(path), "results", "correlated", "binary_snr")
tab_path <- if (length(args) >= 2) args[2] else
  file.path(dirname(path), "results", "thesis_tables")
dir.create(tab_path, showWarnings = FALSE, recursive = TRUE)

population <- readRDS(file.path(res_path, "bin_snr_population.RDS"))

#' x to `digits` decimals, a dash where NA
num <- function(x, digits) {
  ifelse(is.na(x), "---", formatC(x, format = "f", digits = digits))
}

population <- population %>% arrange(scenario, hte_mult)

tab <- tibble(
  hte = paste0("$\\times ", population$hte_mult, "$"),
  bW = num(population$bW, 3),
  ate = num(population$ate, 3),
  sd_tau = num(population$sd_tau, 3),
  sd_noise = num(population$sd_noise, 3),
  snr = num(population$snr, 3),
  clamped = sprintf("%.1f\\%%", 100 * population$share_clamped)
)
group <- label_factor(population$scenario, SS_SCENARIO_LABELS)

out_file <- file.path(tab_path, "bin_snr_population.tex")
grouped_longtable(
  tab,
  col_names = c("HTE multiplier", "$b_W$", "ATE", "SD($\\tau$)", "SD(noise)",
                "SNR", "Treated risk clamped"),
  group = group,
  block = NULL,
  align = c("l", rep("r", ncol(tab) - 1)),
  caption.short = "Population truth, binary signal-to-noise check",
  caption = paste0(
    "Population truth of the binary signal-to-noise check, covariates ",
    "correlated at $\\rho = 0.5$, per scenario and heterogeneity multiplier. ",
    "Values are from the clamped truth over $10^6$ covariate draws. ",
    "HTE multiplier: the factor applied to the effect-modifier coefficients. ",
    "$b_W$: the treatment-model intercept, calibrated at $n = 500$ in every ",
    "cell, so it does not depend on $n$. ATE and SD($\\tau$): the mean and ",
    "standard deviation of the CATE on the risk-difference scale. SD(noise): ",
    "$\\sqrt{\\mathrm{E}[p(1 - p)]}$, the Bernoulli noise averaged over the two ",
    "arms. SNR: SD($\\tau$) / SD(noise). Treated risk clamped: the share of ",
    "the population whose treated risk is held to $[0.01, 0.99]$, with the ",
    "truth clamped to match. The ATE departs from $-0.118$ only where the ",
    "clamp binds."
  ),
  label = "bin_snr_population",
  landscape = FALSE
) %>%
  writeLines(out_file)
message("wrote ", out_file)

# ---- competing risks -----------------------------------------------------------
# competing_risk/snr/ saves no population file: its design (hazard ratios and
# population SNR at rho = 0.5, from rel_eff_snr.R's definition) is the table in
# competing_risk/snr/README.md, copied here. Update both together.

surv <- tibble(
  hte_mult = c(1, 1.5, 2, 3),
  e1_hr = c("1.42 / 2.86", "1.42 / 4.06", "1.42 / 5.75", "1.42 / 11.6"),
  e2_hr = c("0.70 / 0.35", "0.70 / 0.25", "0.70 / 0.17", "0.70 / 0.09"),
  snr1 = c(0.29, 0.43, 0.57, 0.79),
  snr2 = c(0.13, 0.17, 0.20, 0.23)
)

tab <- tibble(
  hte = paste0("$\\times ", surv$hte_mult, "$"),
  e1_hr = surv$e1_hr,
  e2_hr = surv$e2_hr,
  snr1 = num(surv$snr1, 2),
  snr2 = num(surv$snr2, 2)
)

out_file <- file.path(tab_path, "surv_snr_population.tex")
grouped_longtable(
  tab,
  col_names = c("HTE multiplier", "Event 1", "Event 2", "Event 1", "Event 2"),
  group = factor(rep("Scenario 7", nrow(tab))),
  block = NULL,
  header_above = c(" " = 1, "Hazard ratio, $X_3 = 0 / 1$" = 2,
                   "SNR" = 2),
  align = c("l", rep("r", ncol(tab) - 1)),
  caption.short = "Population truth, competing-risk signal-to-noise check",
  caption = paste0(
    "Population truth of the competing-risk signal-to-noise check (scenario ",
    "7, covariates correlated at $\\rho = 0.5$, censoring on), per ",
    "heterogeneity multiplier. HTE multiplier: the factor applied to the ",
    "$X_3$ effect modifiers of both events' hazards; the treatment effect ",
    "terms are unchanged. Hazard ratio: the treatment's cause-specific hazard ",
    "ratio at $X_3 = 0$ and $X_3 = 1$. SNR: SD($\\tau$) / SD(noise) for the ",
    "restricted mean time lost to each event. SNR is a population quantity, ",
    "so it does not change with $n$."
  ),
  label = "surv_snr_population",
  landscape = FALSE
) %>%
  writeLines(out_file)
message("wrote ", out_file)
