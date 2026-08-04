##########
# title: why did the profiling sampler fail to start?
##########
# Throwaway diagnostic. cf_profile.sh subjobs all died with:
#
#   Error in `rs_init(self, private, super, options, wait, wait_timeout)`:
#   ! Could not start R session, timed out
#   1. syrup::syrup({ ...
#   2. callr::r_session$new()
#
# syrup() calls callr::r_session$new() with no arguments, so wait_timeout takes
# its default of 3000 ms, and syrup exposes no way to raise it (its only
# parameters are expr, interval, peak, env). The hypothesis is that spawning a
# fresh R session off /rds/general/project - a networked filesystem - takes
# longer than that, especially once the study's package stack is loaded.
#
# This times R subprocess startup four ways to find out where the cost is.
# Run it through jobscripts/cf_diagnose.sh so the environment matches the
# failure exactly. Takes under a minute.

library(here)

for (pkg in c("callr", "processx", "parallelly", "future")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("diagnostic needs the ", pkg, " package", call. = FALSE)
  }
}

cat("=== environment ===\n")
cat(R.version.string, "\n")
cat("R.home:          ", R.home(), "\n")
cat("libPaths:        ", paste(.libPaths(), collapse = "\n                  "), "\n")
cat("working dir:     ", getwd(), "\n")
cat("availableCores:  ", parallelly::availableCores(), "\n")
cat("callr version:   ", as.character(utils::packageVersion("callr")), "\n")
cat("syrup installed: ", requireNamespace("syrup", quietly = TRUE), "\n")
cat("callr default wait_timeout: 3000 ms\n")

# ---------------------------------------------------------------------------
# the failing call, timed. anything at or above 3s is the explanation.
time_callr <- function(label, k = 3) {
  cat(sprintf("\n--- %s ---\n", label))
  times <- numeric(k)
  for (j in seq_len(k)) {
    t <- system.time(s <- try(callr::r_session$new(), silent = TRUE))
    times[j] <- t[["elapsed"]]
    failed <- inherits(s, "try-error")
    cat(sprintf("  attempt %d: %6.2fs  %s\n", j, times[j],
                if (failed) "FAILED" else "ok"))
    if (!failed) try(s$close(), silent = TRUE)
  }
  cat(sprintf("  median %.2fs, max %.2fs %s\n", median(times), max(times),
              if (max(times) >= 3) "  <-- at or over the 3s budget" else ""))
  invisible(times)
}

cat("\n\n=== 1. shell-level baseline (no R-level machinery) ===\n")
baseline <- system.time(
  system2(file.path(R.home("bin"), "Rscript"), c("-e", "invisible(1)"),
          stdout = NULL, stderr = NULL)
)
cat(sprintf("  Rscript -e 'invisible(1)': %.2fs\n", baseline[["elapsed"]]))

cat("\n\n=== 2. callr before the study stack is loaded ===\n")
cold <- time_callr("cold")

cat("\n\n=== 3. callr after loading the study stack ===\n")
cat("loading cf_models.R (dplyr, future, SuperLearner, gam, coin, GenericML, mlr3, ...)\n")
load_time <- system.time(source(here("crossfitting", "cf_models.R")))
cat(sprintf("  stack load took %.2fs\n", load_time[["elapsed"]]))
warm <- time_callr("warm")

cat("\n\n=== 4. callr with 4 multisession workers running (the failing condition) ===\n")
future::plan(future::multisession, workers = 4)
with_workers <- time_callr("with workers")
future::plan(future::sequential)

cat("\n\n=== 5. the alternatives, for comparison ===\n")
psock <- system.time({
  cl <- parallelly::makeClusterPSOCK(1)
  parallel::stopCluster(cl)
})
cat(sprintf("  parallelly::makeClusterPSOCK(1) [120s timeout]: %.2fs\n",
            psock[["elapsed"]]))

spawn <- system.time({
  p <- processx::process$new(file.path(R.home("bin"), "Rscript"),
                             c("-e", "Sys.sleep(30)"))
})
cat(sprintf("  processx::process$new() [no readiness handshake]: %.2fs\n",
            spawn[["elapsed"]]))
# how long until that process is actually usable, as opposed to merely spawned
alive_at <- Sys.time()
while (!p$is_alive() && as.numeric(Sys.time() - alive_at, units = "secs") < 30) {
  Sys.sleep(0.05)
}
cat(sprintf("  processx process alive: %s\n", p$is_alive()))
p$kill()

# ---------------------------------------------------------------------------
cat("\n\n=== reading ===\n")

worst_callr <- max(c(cold, warm, with_workers))

if (baseline[["elapsed"]] >= 3) {
  cat("The shell-level Rscript baseline alone exceeds 3s. This node cannot spawn\n")
  cat("R processes at a usable rate. FIX: fork-based sampler via parallel::mcparallel\n")
  cat("(no exec, so startup is near-instant) - Linux only, no local Windows testing.\n")
} else if (worst_callr < 1) {
  cat("callr started quickly in every condition. The timeout was not systematic -\n")
  cat("likely a transient node problem. FIX: re-run one profiling index as-is; if it\n")
  cat("recurs, wrap the syrup call in a 3x retry rather than rewriting the sampler.\n")
} else if (max(cold) < 3 && worst_callr >= 3) {
  cat("callr is fast cold but slow once the stack is loaded and/or workers are\n")
  cat("running - memory pressure and contention, not the filesystem alone.\n")
  cat("FIX: processx sampler (crossfitting/cf_sampler.R) - no readiness handshake,\n")
  cat("so nothing can time out.\n")
} else {
  cat("callr exceeds the 3s budget even cold - networked filesystem cost.\n")
  cat("FIX: processx sampler (crossfitting/cf_sampler.R) - no readiness handshake,\n")
  cat("so nothing can time out.\n")
}

cat(sprintf("\nworst callr startup seen: %.2fs against a hardcoded 3.00s budget\n",
            worst_callr))
cat("\ndiagnostic complete\n")
