##########
# title: the case-study learners, scored by the DR score - model evaluation study
##########
# Supplementary analysis (ADEMP.md, "The case-study learners"): does the DR
# score order the applied chapter's three learners - DR-RF (drrf), the causal
# forest (cf) and DR-SL (SL1) - the way their true PEHE does, in this study's
# simulated trials?
#
# Nothing is fit, and no tree is written to:
#   - data, truth, SL1 and the xgb x cv_shared scoring nuisance come from the
#     strategies tree;
#   - drrf and cf are the sample-size study's production fits of the same
#     dataset (me_config.R, CS_LEARNERS), checked to be the same dataset run by
#     run below.
# All three are scored by me_per_model() against the one nuisance, so they
# share phi/phi05 and SL1's scores must equal the strategies tree's exactly.
#
#   Rscript me_cs_learners.R
#
# Writes me_cs_metrics.RDS to CS_RES_PATH: one row per (scenario, n, run,
# model) with true_pehe, mu1_mse, dr and dr05. The thesis tables come from
# results_processing/thesis_tables/me_cs_tables.R. Exits non-zero, after
# saving, if any check below fails.

library(dplyr)
library(here)

source(here("R", "utils.R"))
source(here("model_evaluation", "me_config.R"))
source(here("model_evaluation", "me_metrics.R")) # me_per_model()

# the strategies tree's names for the two scores kept here
score_col <- function(type) paste(type, CS_ARM, CS_PIPELINE, sep = "_")
SCORE_COLS <- c(dr = score_col("dr"), dr05 = score_col("dr05"))

#' One simulated trial's three learners, scored; NULL for a run the main study
#' excluded (no strategies file)
#'
#' @param i a row of study_strat$grid, and so of study_ss_cts$grid, which is
#'   the same grid with rho added
score_run <- function(i) {
  param <- study_strat$grid[i, ]
  f <- paste0("res_sim_", param$run, ".RDS")
  strat_f <- file.path(combo_dir(study_strat, param), f)
  if (!file.exists(strat_f)) return(NULL)
  ss_f <- file.path(combo_dir(study_ss_cts, study_ss_cts$grid[i, ]), f)
  if (!file.exists(ss_f)) stop("no sample-size run at ", ss_f)

  strat <- readRDS(strat_f)
  ss <- readRDS(ss_f)

  # The same simulated trial. This study's X is the sample-size study's
  # standardised (prepare_design_matrix()), so X is compared after scale().
  tag <- sprintf("scenario %s, n %s, run %s", param$scenario, param$n, param$run)
  same <- c(
    Y = identical(strat$data$Y, ss$data$Y),
    W = identical(strat$data$W, ss$data$W),
    X = isTRUE(all.equal(unname(as.matrix(strat$data$X)),
                         unname(scale(as.matrix(ss$data[, -(1:2)]))),
                         check.attributes = FALSE)),
    tau = identical(strat$truth$tau, ss$truth$tau)
  )
  if (!all(same)) {
    stop(tag, ": not the sample-size study's dataset (",
         paste(names(same)[!same], collapse = ", "), " differ)")
  }

  learners <- lapply(CS_LEARNERS, function(field) {
    fit <- if (field == "SL1") strat$SL1 else ss[[field]]
    if (is.null(fit$tau)) stop(tag, ": no ", field, "$tau")
    list(tau = fit$tau)
  })

  nuis <- strat$nuisances[[CS_PIPELINE]][[CS_ARM]]
  if (is.null(nuis)) stop(tag, ": no ", CS_PIPELINE, " ", CS_ARM, " nuisance")

  sim_res <- c(learners, list(
    data = strat$data,
    truth = strat$truth,
    nuisances = setNames(list(setNames(list(nuis), CS_ARM)), CS_PIPELINE)
  ))
  keys <- tibble(scenario = param$scenario, n = param$n, run = param$run)

  bind_rows(lapply(names(CS_LEARNERS), function(m) {
    s <- me_per_model(sim_res[[m]], sim_res$truth$tau, m, sim_res, keys)
    bind_cols(keys, tibble(
      model = m, true_pehe = s$true_pehe, mu1_mse = s$mu1_mse,
      dr = s[[SCORE_COLS[["dr"]]]], dr05 = s[[SCORE_COLS[["dr05"]]]]
    ))
  }))
}

grid <- study_strat$grid
rows <- lapply(seq_len(nrow(grid)), function(i) {
  if (i %% 100 == 0) message(i, "/", nrow(grid))
  score_run(i)
})
metrics <- bind_rows(rows)
if (nrow(metrics) == 0) {
  stop("no strategies-tree runs found under ", study_strat$res_path)
}

dir.create(CS_RES_PATH, showWarnings = FALSE, recursive = TRUE)
saveRDS(metrics, file.path(CS_RES_PATH, "me_cs_metrics.RDS"))
message("wrote ", file.path(CS_RES_PATH, "me_cs_metrics.RDS"))

# ---- checks ------------------------------------------------------------------

fails <- character()
key_cols <- c("scenario", "n", "run", "model")

# completeness, as me_results.qmd reports it: runs per cell, and any run
# without all three learners scored
excluded <- grid[vapply(rows, is.null, logical(1)), c("scenario", "n", "run")]
cat(sprintf("\n%d simulated trials scored; %d excluded from the main study:\n",
            nrow(distinct(metrics, scenario, n, run)), nrow(excluded)))
if (nrow(excluded)) print(as.data.frame(excluded), row.names = FALSE)
cat("\nsimulated trials per cell:\n")
print(as.data.frame(count(distinct(metrics, scenario, n, run), scenario, n,
                          name = "trials")),
      row.names = FALSE)

incomplete <- metrics %>%
  group_by(scenario, n, run) %>%
  filter(n() != length(CS_LEARNERS) |
           any(!is.finite(c(true_pehe, dr, dr05)))) %>%
  ungroup()
if (nrow(incomplete)) {
  fails <- c(fails, sprintf("%d row(s) in runs with a missing or non-finite score",
                            nrow(incomplete)))
}

# SL1 must come back exactly as the strategies tree scored it, and its true
# PEHE as the main tree did
sl1 <- filter(metrics, model == "SL1")
check_sl1 <- function(file, cols) {
  if (!file.exists(file)) {
    cat("\nskipped the SL1 check against", file, "(not found)\n")
    return(invisible())
  }
  ref <- readRDS(file) %>%
    filter(model == "SL1") %>%
    select(all_of(key_cols), all_of(unname(cols)))
  names(ref)[match(cols, names(ref))] <- paste0("ref_", names(cols))
  joined <- inner_join(sl1, ref, by = key_cols)
  if (nrow(joined) != nrow(sl1)) {
    fails <<- c(fails, sprintf("SL1: %d of %d runs found in %s",
                               nrow(joined), nrow(sl1), basename(file)))
  }
  for (nm in names(cols)) {
    same <- identical(joined[[nm]], joined[[paste0("ref_", nm)]])
    cat(sprintf("SL1 %-9s identical to %s: %s\n", nm, basename(file), same))
    if (!same) fails <<- c(fails, sprintf("SL1 %s differs from %s", nm, basename(file)))
  }
}
cat("\n")
check_sl1(file.path(study_strat$res_path, "me_strat_metrics.RDS"),
          c(true_pehe = "true_pehe", mu1_mse = "mu1_mse",
            dr = SCORE_COLS[["dr"]], dr05 = SCORE_COLS[["dr05"]]))
# the main tree has only the `whole` arm, so only the truth-based columns
check_sl1(file.path(study$res_path, "me_metrics.RDS"),
          c(true_pehe = "true_pehe", mu1_mse = "mu1_mse"))

# drrf and cf against the sample-size study's mean MSE over its 500 runs
# (continuous, rho = 0.5). Close, not equal: these are its runs 1-100.
ss_ref <- tibble::tribble(
  ~scenario, ~n,   ~drrf, ~cf,
  1,         250,  0.018, 0.008,
  1,         500,  0.011, 0.005,
  1,         1000, 0.008, 0.003,
  2,         250,  0.135, 0.187,
  2,         500,  0.070, 0.089,
  2,         1000, 0.039, 0.045,
  3,         250,  0.291, 0.447,
  3,         500,  0.200, 0.276,
  3,         1000, 0.129, 0.178,
  4,         250,  0.101, 0.147,
  4,         500,  0.056, 0.072,
  4,         1000, 0.030, 0.034
)
cat("\nmean true PEHE here (runs 1-100) vs the sample-size study's MSE (500 runs):\n")
here_means <- metrics %>%
  filter(model %in% c("drrf", "cf")) %>%
  group_by(scenario, n, model) %>%
  summarise(mean = mean(true_pehe), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = model, values_from = mean)
print(as.data.frame(
  inner_join(here_means, ss_ref, by = c("scenario", "n"), suffix = c("", "_ss")) %>%
    transmute(scenario, n,
              drrf = round(drrf, 3), drrf_ss,
              cf = round(cf, 3), cf_ss) %>%
    arrange(scenario, n)
), row.names = FALSE)

if (length(fails)) {
  cat("\nFAILURES:\n")
  cat(paste0("  - ", fails, collapse = "\n"), "\n")
  quit(save = "no", status = 1)
}
cat("\nAll checks passed.\n")
