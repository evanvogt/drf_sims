##########
# title: application sizing - the one definition of the timing grid
##########
# Sourced by as_analysis.R, as_check.R, as_collect.R and as_summary.R. The array
# index is a row number of `grid`, so `grid` must never be filtered or
# reordered; add rows at the end.
#
# One row is one timing job. `task` says what it times:
#   impute  handle_missingness() on a whole platform domain (all arms, by arm)
#           or on the two-arm dataset, for one handling method. Nothing is
#           analysed; the analysis jobs below time a completed dataset instead.
#   binary  platform comparison: the 6 binary outcomes (bin1, bin2, and
#           discharge by each horizon), one cate_methods() fit each
#   surv    platform comparison: all_cate_surv_models() at each of the 4
#           horizons, RMTL for discharge only
#   cts     two-arm: the 3 continuous outcomes, one cate_methods() fit each
#   ci      one outcome, cate_methods() with the half-sample bootstrap on top
#   sf      one outcome, find_optimal_sf() at ONE sample.fraction (`sf`); a
#           full calibration is the sum over the 10 values
#
# `handling` on the analysis tasks is the shape of the analysed data, not the
# imputation that produced it - an analysis takes as long on a mean-imputed
# dataset as on a missForest one:
#   completed       every covariate filled in (mean imputation stands in for
#                   all of mean / missForest / regression / each of the MI
#                   datasets)
#   indicator       filled in plus one indicator per incomplete covariate
#   none            NAs left in (grf arms only; SuperLearner arms and the RSF
#                   arms cannot take them)
#   complete_cases  two-arm only: the complete rows
#   IPW             two-arm only: the complete rows, weighted
# so as_summary.R builds each method's total from its impute row and the
# matching analysis rows (multiple imputation = impute + 50 x completed).

library(here)
source(here("R", "pipeline.R"))
source(here("application_sizing", "as_schema.R"))

# ---- analysis settings ------------------------------------------------------

N_FOLDS <- 10
N_IMP <- 50

# no oracle-type arms
CATE_MODELS_APPLIED <- c("causal_forest", "dr_random_forest", "dr_superlearner")
SURV_MODELS_APPLIED <- c(
  "csf_sh",
  "pseudo_cf_whole_oob", "pseudo_cf_whole_scf", "pseudo_cf_cvps_scf",
  "pseudo_dr_whole_oob", "pseudo_dr_whole_scf", "pseudo_dr_cvps_scf",
  "sl_dr_whole", "sl_dr_cvps",
  "rsf_dr_oob", "rsf_dr_scf"
)
# the ones that run with NAs in X (grf handles them; SuperLearner and rfsrc
# do not)
SURV_MODELS_NA <- c(
  "csf_sh",
  "pseudo_cf_whole_oob", "pseudo_cf_whole_scf", "pseudo_cf_cvps_scf",
  "pseudo_dr_whole_oob", "pseudo_dr_whole_scf", "pseudo_dr_cvps_scf"
)
SURV_ESTIMANDS_APPLIED <- "RMTL1"   # discharge

# Outcomes used as predictors when imputing a domain / the dataset. The second
# binary outcome is left out of the platform set: it is the death indicator,
# which `status` already carries, and an exact copy adds nothing but a
# rank-deficient regression imputation.
IMPUTE_OUTCOMES <- list(platform = c("bin1", "time", "status"),
                        two_arm  = c("y1", "y2", "y3"))

# bootstrap CIs: the CI studies' settings; sf at the top of the calibration
# grid, the slowest
CI <- list(boot = 200, sf = 0.5, alpha = 0.05)
SF_GRID <- seq(0.05, 0.5, 0.05)
SF_N_SIM <- 50
SF_CI_BOOT <- 100
CI_OUTCOME <- c(platform = "bin1", two_arm = "y1")

PLATFORM_HANDLING <- c("mean_imputation", "missforest", "regression",
                       "missing_indicator", "multiple_imputation")
TWO_ARM_HANDLING <- c("complete_cases", "IPW", PLATFORM_HANDLING)

# ---- the grid ---------------------------------------------------------------

row_block <- function(dataset, task, handling, sf = 0) {
  expand.grid(dataset = dataset, task = task, handling = handling, sf = sf,
              run = 1, stringsAsFactors = FALSE)
}

cmp_ids <- platform_comparisons()$id

grid <- rbind(
  row_block(names(PLATFORM_DOMAINS), "impute", PLATFORM_HANDLING),
  row_block("two_arm", "impute", TWO_ARM_HANDLING),
  row_block(cmp_ids, c("binary", "surv"), c("completed", "indicator", "none")),
  row_block("two_arm", "cts", c("completed", "indicator", "none",
                                "complete_cases", "IPW")),
  row_block(c(cmp_ids, "two_arm"), "ci", "completed"),
  row_block(c(cmp_ids, "two_arm"), "sf", "completed", SF_GRID)
)
rownames(grid) <- NULL

study <- study_config(
  name     = "application_sizing",
  prefix   = "as",
  res_path = file.path(dirname(here()), "results", "application_sizing"),
  grid     = grid,
  path_cols   = c("dataset", "task", "handling", "sf"),
  path_prefix = c(sf = "sf_"),
  n_sims      = 1,
  failed_file = here("application_sizing", "jobscripts", "failed_ids.txt")
)
