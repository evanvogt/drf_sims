##########
# title: application sizing - what the full analysis would cost
##########
# Rscript as_summary.R [speed_ratio]
#
# Reads as_collect.R's CSVs and builds, per trial and missing-data method, the
# total run time of the whole analysis, plus the cost of CIs and of an sf
# calibration per outcome. `speed_ratio` (default 1) scales cluster seconds to
# the laptop: time one cheap row on both and pass laptop / cluster.
#
# How a method's total is assembled from the rows (see as_config.R):
#   mean / missForest / regression   impute + analysis on "completed" data
#   missing indicator                impute + analysis on "indicator" data
#   none                             analysis on "none" data
#   multiple imputation              impute + N_IMP x analysis on "completed"
#   complete cases / IPW (two-arm)   impute + their own analysis rows

library(here)
source(here("application_sizing", "as_config.R"))

args <- as.numeric(commandArgs(trailingOnly = TRUE))
speed <- if (length(args) >= 1 && !is.na(args[1])) args[1] else 1

fits <- read.csv(file.path(study$res_path, "as_fits.csv"))
jobs <- read.csv(file.path(study$res_path, "as_jobs.csv"))
fits$hours <- fits$seconds * speed / 3600
fits$trial <- ifelse(fits$dataset == "two_arm", "two_arm", "platform")

hrs <- function(trial, task, handling = NULL) {
  keep <- fits$trial == trial & fits$task %in% task
  if (!is.null(handling)) keep <- keep & fits$handling %in% handling
  sum(fits$hours[keep])
}

shape_of <- c(mean_imputation = "completed", missforest = "completed",
              regression = "completed", missing_indicator = "indicator",
              multiple_imputation = "completed", none = "none",
              complete_cases = "complete_cases", IPW = "IPW")

method_table <- function(trial) {
  methods <- c(if (trial == "platform") PLATFORM_HANDLING else TWO_ARM_HANDLING, "none")
  tasks <- if (trial == "platform") c("binary", "surv") else "cts"
  do.call(rbind, lapply(methods, function(m) {
    imp <- if (m == "none") 0 else hrs(trial, "impute", m)
    analysis <- hrs(trial, tasks, shape_of[[m]])
    reps <- if (m == "multiple_imputation") N_IMP else 1
    data.frame(trial = trial, method = m, impute_h = imp,
               analysis_h = analysis * reps, total_h = imp + analysis * reps)
  }))
}

cat("Hours on", if (speed == 1) "the cluster" else paste0("the laptop (x", speed, ")"), "\n\n")

totals <- rbind(method_table("platform"), method_table("two_arm"))
print(format(totals, digits = 3), row.names = FALSE)

# CIs: the fit with the bootstrap replaces the plain fit, so the extra cost per
# outcome is the difference (summed over the platform's comparisons)
plain <- function(trial) {
  out <- CI_OUTCOME[[trial]]
  task <- if (trial == "platform") "binary" else "cts"
  sum(fits$hours[fits$trial == trial & fits$task == task &
                   fits$handling == "completed" & fits$fit == out])
}
ci <- data.frame(
  trial = c("platform", "two_arm"),
  ci_fit_h = c(hrs("platform", "ci"), hrs("two_arm", "ci")),
  plain_fit_h = c(plain("platform"), plain("two_arm")),
  sf_calibration_h = c(hrs("platform", "sf"), hrs("two_arm", "sf"))
)
ci$ci_extra_per_outcome_h <- ci$ci_fit_h - ci$plain_fit_h
cat("\nPer outcome (platform: summed over its", length(platform_comparisons()$id),
    "comparisons):\n")
print(format(ci, digits = 3), row.names = FALSE)

cat("\nPeak memory by task (GB, R-side upper bound; PBS peak where available):\n")
mem <- aggregate(cbind(peak_r_gb = peak_r_mb / 1024) ~ task, data = jobs, FUN = max)
if ("pbs_peak_gb" %in% names(jobs) && any(!is.na(jobs$pbs_peak_gb))) {
  mem <- merge(mem, aggregate(pbs_peak_gb ~ task, data = jobs, FUN = max), all.x = TRUE)
}
print(format(mem, digits = 3), row.names = FALSE)

write.csv(totals, file.path(study$res_path, "as_method_totals.csv"), row.names = FALSE)
write.csv(ci, file.path(study$res_path, "as_ci_costs.csv"), row.names = FALSE)
