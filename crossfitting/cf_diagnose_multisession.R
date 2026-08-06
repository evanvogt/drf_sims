##########
# title: why don't multisession workers show up in the profiling process tree?
##########
# Companion to cf_diagnose_sampler.R (same failure family: R subprocess startup
# off /rds/general/project - a networked filesystem - being slow enough to miss a
# default timeout). That diagnostic found callr::r_session$new() (syrup's sampler)
# timing out at its hardcoded 3s budget. This checks whether future::multisession's
# PSOCK workers hit the same wall - which would explain cf_profile.R's
# "expected N processes... saw 1" warning on workers > 1 cells: not
# oversubscription, not a syrup ps-detection blind spot, but the workers never
# actually connecting.
#
# Run through jobscripts/cf_diagnose.sh so the environment matches cf_profile.sh
# exactly. Three conditions, each: establish a plan, dispatch one future per
# worker, and independently verify via ps::ps_children() (not future's own
# bookkeeping, which would just report what it believes happened) whether that
# many child processes actually exist.

library(here)
library(future)
library(furrr)

for (pkg in c("parallelly", "ps")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("diagnostic needs the ", pkg, " package", call. = FALSE)
  }
}

cat("=== environment ===\n")
cat(R.version.string, "\n")
cat("availableCores:     ", parallelly::availableCores(), "\n")
cat("supportsMulticore:  ", parallelly::supportsMulticore(), "\n")
cat("future version:     ", as.character(utils::packageVersion("future")), "\n")
cat("parallelly version: ", as.character(utils::packageVersion("parallelly")), "\n")

workers <- 4

# ps_children() returns a list of ps_handle objects, not a data.frame - nrow()
# on that silently gives NULL rather than erroring, so use length() instead.
count_children <- function() length(ps::ps_children(ps::ps_handle(), recursive = TRUE))

run_condition <- function(label, plan_call) {
  cat(sprintf("\n--- %s ---\n", label))
  t <- system.time({
    plan_call()
    pids <- future_map(seq_len(workers), function(i) Sys.getpid(),
                       .options = furrr_options(seed = TRUE))
  })
  n_seen <- count_children()
  n_distinct_worker_pids <- length(unique(unlist(pids)))
  plan(sequential)
  cat(sprintf("  setup+dispatch: %.2fs\n", t[["elapsed"]]))
  cat(sprintf("  distinct worker PIDs future reports: %d (expected %d)\n",
              n_distinct_worker_pids, workers))
  cat(sprintf("  child processes ps sees post-dispatch: %d (expected %d)\n",
              n_seen, workers))
  list(label = label, elapsed = t[["elapsed"]],
       future_pids = n_distinct_worker_pids, ps_children = n_seen)
}

results <- list()

results$default <- run_condition("1. multisession, defaults", function() {
  plan(multisession, workers = workers)
})

results$raised_timeout <- run_condition("2. multisession, raised PSOCK timeouts", function() {
  old <- options(parallelly.makeNodePSOCK.connectTimeout = 120,
                 parallelly.makeNodePSOCK.timeout = 120)
  on.exit(options(old), add = TRUE)
  plan(multisession, workers = workers)
})

if (parallelly::supportsMulticore()) {
  results$multicore <- run_condition("3. multicore (fork)", function() {
    plan(multicore, workers = workers)
  })
} else {
  cat("\n--- 3. multicore (fork) ---\n  skipped: not supported on this platform\n")
}

cat("\n\n=== reading ===\n")
if (results$default$ps_children == workers) {
  cat("multisession defaults already produce the expected worker count here - the\n")
  cat("original warning may have been a one-off or specific to a different node/cell.\n")
  cat("FIX: re-run one profiling index as-is; if it recurs, escalate rather than\n")
  cat("changing plan() calls blind.\n")
} else if (results$raised_timeout$ps_children == workers) {
  cat("raising parallelly's PSOCK timeout fixes it - the workers were timing out\n")
  cat("during connection setup, same failure family as the callr/syrup issue.\n")
  cat("FIX: set parallelly.makeNodePSOCK.connectTimeout/.timeout before plan()\n")
  cat("calls in cf_profile.R, cf_analysis.R, cf_testing.R.\n")
} else if (!is.null(results$multicore) && results$multicore$ps_children == workers) {
  cat("raising the timeout did NOT fix it, but multicore (fork) does.\n")
  cat("FIX: switch plan(multisession, ...) to plan(multicore, ...) - guarded by\n")
  cat("parallelly::supportsMulticore() - in cf_profile.R, cf_analysis.R,\n")
  cat("cf_testing.R.\n")
} else {
  cat("none of the three conditions produced the expected worker count. this isn't\n")
  cat("a startup-timeout problem - something else is blocking PSOCK/fork workers\n")
  cat("on this node. escalate rather than guessing further.\n")
}

cat("\ndiagnostic complete\n")
