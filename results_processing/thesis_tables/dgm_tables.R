##########
# title: appendix tables - the true CATE structure of the DGMs
##########
# The numbers behind thesis_figures/dgm.R, from the same draws
# (results_processing/dgm_population.R). No simulation output is read.
#
# Writes to ../results/thesis_tables/:
#   dgm_cts_bin.tex   continuous_corr_0.5 and binary_corr_0.5, scenarios 1-4:
#                     SD(tau), the central 98% of tau - ATE, cor(m0, tau), the
#                     R^2 of a linear fit of tau, and at each n the ATE with
#                     the share of the population whose CATE has the opposite
#                     sign
#   dgm_surv.tex      competing risks, rho = 0.5, scenario x event: the RMTL
#                     CATE's mean and SD, its mean by X3, and X3's effect with
#                     X1 and X2 held (surv_dgm_check.R's standardised
#                     separation)
#
# Local, a few minutes (the competing-risk truth is numerical integration):
#   cd results_processing/thesis_tables; Rscript dgm_tables.R

library(here)
source(here("R", "figures.R"))
source(here("R", "tables.R"))
source(here("results_processing", "dgm_population.R"))

# paths
path <- here()
tab_path <- file.path(dirname(path), "results", "thesis_tables")
dir.create(tab_path, showWarnings = FALSE, recursive = TRUE)

OUTCOME_LABELS <- c(continuous = "Continuous (mean difference)",
                    binary = "Binary (risk difference)")

write_table <- function(tex, out) {
  out_file <- file.path(tab_path, paste0(out, ".tex"))
  writeLines(tex, out_file)
  message("wrote ", out_file)
}

#' x to `digits` decimals (one per element, or one for all), a dash where NA
num <- function(x, digits) {
  out <- mapply(function(v, d) formatC(v, format = "f", digits = d), x, digits)
  ifelse(is.na(x), "---", out)
}

# ---- continuous / binary ------------------------------------------------------

cv <- dgm_covariates()
summ <- cate_summary(cate_population(cv), cv)

# risk differences need a third decimal; mean differences do not
dig <- ifelse(summ$outcome == "binary", 3, 2)
cts_bin <- tibble(
  scenario = as.character(label_factor(summ$scenario, SS_SCENARIO_LABELS)),
  sd = num(summ$sd_tau, dig),
  range = ifelse(summ$sd_tau > 0,
                 paste0("$[", num(summ$q01, dig), ", ", num(summ$q99, dig), "]$"),
                 "---"),
  cor = num(summ$cor_m0_tau, 2),
  r2 = num(summ$r2_linear, 2),
  !!!setNames(lapply(DGM_N_LEVELS, function(n) {
    paste0(num(summ[[paste0("ate_", n)]], dig), " (",
           sprintf("%.0f", 100 * summ[[paste0("opp_", n)]]), "\\%)")
  }), paste0("n", DGM_N_LEVELS))
)
group <- label_factor(summ$outcome, OUTCOME_LABELS)

grouped_longtable(
  cts_bin,
  col_names = c("Scenario", "SD($\\tau$)", "$\\tau - $ ATE, 1--99\\%",
                "cor($m_0$, $\\tau$)", "Linear $R^2$",
                paste0("$n = ", DGM_N_LEVELS, "$")),
  group = group,
  block = NULL,
  # no % here: add_header_above() unescapes \%, which would comment out the line
  header_above = c(" " = 5, "ATE (share with the opposite sign)" = length(DGM_N_LEVELS)),
  align = c("l", rep("r", ncol(cts_bin) - 1)),
  caption.short = "True CATE structure, continuous and binary outcomes",
  caption = paste0(
    "True CATE structure of the continuous and binary data-generating ",
    "mechanisms, scenarios 1--4, covariates correlated at $\\rho = 0.5$. ",
    "SD($\\tau$) is exact, by quadrature; the other columns are over ",
    format(DGM_N_POP, big.mark = ",", scientific = FALSE), " covariate draws. ",
    "$\\tau - $ ATE, 1--99\\%: the 1st and 99th percentiles of the CATE about ",
    "the ATE, which do not depend on $n$. cor($m_0$, $\\tau$): the correlation ",
    "of the control mean (binary: control risk) with the CATE. Linear $R^2$: ",
    "the share of the CATE's variance a linear model in $X_1$--$X_5$ and ",
    "$X_{01}$--$X_{03}$ captures. The ATE is calibrated to 80\\% power at each ",
    "$n$, so it shrinks as $n$ grows while the heterogeneity about it does ",
    "not; in brackets, the share of the population whose CATE has the opposite ",
    "sign to the ATE. The binary effect modifiers are the continuous ones with ",
    "the sign reversed, scaled onto the risk-difference scale. --- : not ",
    "defined (constant CATE)."
  ),
  label = "dgm_cts_bin",
  landscape = FALSE
) %>%
  write_table("dgm_cts_bin")

# ---- competing risks ------------------------------------------------------------

ssumm <- surv_summary(surv_population(cv))

surv <- ssumm %>%
  mutate(
    scenario = label_factor(scenario, SURV_SCENARIO_LABELS),
    event = as.character(event),
    across(c(mean_tau, sd_tau, mean_x3_0, mean_x3_1, x3_effect), ~ num(.x, 2))
  ) %>%
  arrange(scenario, event)

grouped_longtable(
  select(surv, event, mean_tau, sd_tau, mean_x3_0, mean_x3_1, x3_effect),
  col_names = c("Event", "Mean", "SD", "$X_3 = 0$", "$X_3 = 1$",
                "$X_3$ effect"),
  group = surv$scenario,
  block = NULL,
  header_above = c(" " = 3, "Mean, by $X_3$" = 2, " " = 1),
  align = c("l", rep("r", 5)),
  caption.short = "True CATE structure, competing risks",
  caption = paste0(
    "True CATE structure of the competing-risks data-generating mechanism, ",
    "covariates correlated at $\\rho = 0.5$: the CATE on each event's ",
    "restricted mean time lost (RMTL) to the 28-day horizon, in days, over ",
    format(DGM_N_POP_SURV, big.mark = ",", scientific = FALSE), " covariate ",
    "draws. Positive: treatment brings that event forward. The hazards' ",
    "effect modifier is $X_3$ alone; the spread within an $X_3$ group comes from ",
    "the prognostic $X_1$ and $X_2$, through the non-linear RMTL scale and the ",
    "competition between the events. $X_3$ effect: the mean change in the CATE ",
    "when $X_3$ is switched from 0 to 1 with $X_1$ and $X_2$ held, unlike the ",
    "difference of the two by-$X_3$ means, which under correlation also ",
    "picks up the different $X_1$ and $X_2$ of the two groups."
  ),
  label = "dgm_surv",
  landscape = FALSE
) %>%
  write_table("dgm_surv")
