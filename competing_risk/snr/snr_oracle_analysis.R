##########
# Title: competing-risk SNR check - oracle DR-learner on the saved runs
##########
# Reads grid rows <first> .. <first> + <count> - 1 of the main SNR study's
# result files and fits the full-oracle DR-learner on each (snr_oracle_config.R):
#   Rscript snr_oracle_analysis.R <first> [count] [grf_threads]
# count defaults to 1, so a rerun of failed_ids_oracle.txt passes one index.
# jobscripts/snr_oracle_1.sh steps the array by 100, a cell per subjob.
#
# A row whose main result file is missing is skipped with a message and writes
# nothing, so study.R check on snr_oracle_config.R lists it. Rerun the main
# study's failures first. A row that errors does not stop the rest of the
# block; the script exits with status 1 at the end if any did.

library(dplyr)
library(grf)
library(here)

source(here("R", "utils.R"))
source(here("R", "cate_models.R"))
source(here("competing_risk", "surv_dgm.R"))
source(here("competing_risk", "surv_models.R"))
source(here("competing_risk", "snr", "snr_oracle_config.R"))

args <- as.numeric(commandArgs(trailingOnly = T))
first <- args[1]
count <- if (length(args) >= 2 && !is.na(args[2])) args[2] else 1
grf_threads <- if (length(args) >= 3 && !is.na(args[3])) args[3] else 1

horizon <- survival_scenario_params$event_horizon[1]

rows <- first:min(first + count - 1, nrow(study$grid))
failed <- integer(0)

for (i in rows) {
  param <- study$grid[i, ]
  src <- file.path(combo_dir(snr_study, param), paste0("res_sim_", param$run, ".RDS"))
  if (!file.exists(src)) {
    cat(sprintf("row %d: no %s - main run missing, skipped\n", i, src))
    next
  }

  ok <- tryCatch({
    t0 <- Sys.time()
    res <- readRDS(src)
    data <- res$data
    truth <- res$truth
    rm(res)

    # as all_cate_surv_models() builds them
    X <- as.matrix(data[, !(names(data) %in% c("Y", "D", "W"))])
    W <- data$W
    pseudo <- pseudo_all(data$Y, data$D, horizon)
    stopifnot(nrow(truth) == nrow(X))

    # seeds the stage-2 forests, which draw grf's seed from R's RNG
    setup_rng_stream(param$run)
    oracle <- setNames(lapply(SNR_ESTIMANDS, function(e) {
      po <- dr_pseudo(pseudo[[paste0("ps_", e)]], W,
                      truth[[paste0(e, "_1")]], truth[[paste0(e, "_0")]], TRUE_PS)
      stage2_whole_rf(X, po, num.threads = grf_threads)$tau
    }), SNR_ESTIMANDS)

    save_run(setNames(list(oracle), ORACLE_ARM), combo_dir(study, param), param$run)
    cat(sprintf("row %d (scenario %d, n %d, hte_mult %g, run %d): %.1f s\n",
                i, param$scenario, param$n, param$hte_mult, param$run,
                as.numeric(Sys.time() - t0, units = "secs")))
    TRUE
  }, error = function(err) {
    cat(sprintf("row %d FAILED: %s\n", i, conditionMessage(err)))
    FALSE
  })
  if (!ok) failed <- c(failed, i)
}

if (length(failed)) {
  cat("failed rows:", failed, "\n")
  quit(status = 1, save = "no")
}
print("Oracle block completed!")
