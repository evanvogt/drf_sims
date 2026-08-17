##########
# Title: why the array indices in jobscripts/failed_ids.txt fail
##########
# The HPC is only keeping the PBS `.o` files, so the R error that killed these
# runs is invisible on the cluster. This script reproduces it locally. It
# diagnoses only - it changes nothing in surv_models.R.
#
# STATUS: the bug this script found is FIXED (see pseudo_sl_t_split()'s header
# comment in surv_models.R). The script is kept because it is the tool for the
# next batch of missing runs, not because the diagnosis is still open - run it
# against a fresh failed_ids.txt and it will tell you whether you are looking at
# the same two routes or something new. On the fixed code every stage should
# come back clean.
#
# ROOT CAUSE (of the 225 failures of 1400, as found): every failure reproduced
# was in ONE arm, pseudo_sl_t_split(), by TWO different routes. That arm was the
# only SuperLearner call site in this study that passed
# `SL.library = sl_library` unvalidated AND computed its pseudo-values with an
# unguarded pseudoyl()/pseudomean(). Every other SL call goes through
# R/cate_models.R::pretest_superlearner(), and every other pseudo-value path
# carries pseudo_crossfit()'s is.na() -> pseudo_whole fallback. That arm now has
# both guards too.
#
#   MODE A - degenerate library (mostly censoring = FALSE)
#
#   In the scenarios with bW_1 = -0.7 (1, 3, 4, 6, 7 - see
#   survival_scenario_params in surv_dgm.R) almost every treated subject has the
#   cause-1 event well before horizon = 28, so cause-2 events before the horizon
#   collapse in the treated arm: on index 67, 5 of the 500 subjects against 53
#   in the control arm. The training-fold RMTL2 pseudo-values for W == 1 then
#   take only a handful of distinct values (3 across 195 rows, on index 67's
#   fold 5 - stage 2 reports this as min_n_unique). Inside SuperLearner's own
#   5-fold CV at least one fold gets a constant y, so SL.glmnet errors with
#
#     y is constant; gaussian glmnet fails at standardization step
#
#   and SL.glm/SL.gam are dropped for NA/error alongside it. Their $fitLibrary
#   entries are left NULL, and make_t_cate()'s predict(sl1, newdata = X_test)
#   (surv_models.R:836-837) defaults to onlySL = FALSE - so it predicts from
#   EVERY library member, NULL fits included:
#
#     no applicable method for 'predict' applied to an object of class "NULL"
#
#   MODE B - NA pseudo-values (mostly censoring = TRUE)
#
#   The training-fold pseudoyl()/pseudomean() at surv_models.R:818-819 return NA
#   for whoever holds the maximum observed time in that training set - the
#   pseudo::ci.omit bug already written up in README "Known issues". Nothing
#   guards them here, so the NA goes straight into SuperLearner():
#
#     missing data is currently not supported. Check Y, X, and newX
#
#   Confirmed on indices 85, 183, 771 and 897: the fold whose training
#   pseudo-values carry NAs is exactly the fold furrr names in the error.
#   README's "T-learner split pseudo-obs - same latent risk, unverified" called
#   this; it is no longer unverified.
#
#   Either way the run aborts before saveRDS(), which is why check_failed() sees
#   no res_sim_<run>.RDS. Scenarios 2 and 5 (bW_1 = 0) never fail; the scenarios
#   that also carry bW_2 = +0.7 (4, 6, 7) fail most, because that slows the
#   competing event in the treated arm further.
#
#   Stage 2 measures both routes per index and labels which one it predicts, so
#   the split between them is read off the data rather than assumed.
#
# Run from competing_risk/:
#   Rscript surv_failed_diagnose.R                        # stages 1 + 2
#   Rscript surv_failed_diagnose.R --stage=1
#   Rscript surv_failed_diagnose.R --stage=2 --controls=100
#   Rscript surv_failed_diagnose.R --stage=3 --ids=67,68
#   Rscript surv_failed_diagnose.R --stage=4 --logs=jobscripts/logs_1
#   Rscript surv_failed_diagnose.R --stage=4 --logs=jobscripts/logs_rerun --rerun
#
# Costs, measured on R 4.5.3 locally with plan(sequential):
#   stage 1  instant
#   stage 2  ~12 s an index, so ~65 min for the 225 failed plus 100 controls
#            (most of it the N_PROBE_CELLS x SL_ATTEMPTS SuperLearner fits)
#   stage 3  ~15 min an index - the SuperLearner arms are almost all of it -
#            so it defaults to two indices and is meant to be pointed at a
#            handful, not at the whole failed list
# Run stages 2 and 3 backgrounded, e.g.
#   Rscript surv_failed_diagnose.R --stage=2 > stage2.log 2>&1 &
#
# Stage 3's log is interleaved with pretest_superlearner()'s own "Removed
# libraries due to NA/error" cat() - that is left in on purpose. It is the
# guarded arms dropping the same learners, on the same data, that
# pseudo_sl_t_split() goes on to hand the unvalidated library to.
#
# CSVs land in <results>/competing_risk/diagnostics/.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(furrr)
  library(grf)
  library(pseudo)
  library(SuperLearner)
})

# Same source order as surv_analysis.R, so pretest_superlearner, stage2_whole_rf,
# stage_2_sl, dr_pseudo and trim_ps resolve to exactly what the failing runs saw.
source(here("R", "utils.R"))
source(here("R", "cate_models.R"))
source(here("competing_risk", "surv_dgm.R"))
source(here("competing_risk", "surv_models.R"))
source(here("competing_risk", "surv_config.R"))

# production settings (surv_analysis.R)
HORIZON <- 28
ESTIMANDS <- c("RMTL1", "RMTL2", "RMSTc")

# How many of an index's least varied (fold, estimand, arm) cells stage 2 puts
# through a real SuperLearner fit, and how many times each.
#
# Repeats are the point, not belt-and-braces. Whether a degenerate cell trips
# SL.glmnet depends on where SuperLearner's own 5-fold CV happens to cut, and
# this script cannot reproduce the RNG state pseudo_sl_t_split() had. Measured
# over 8 attempts on the least varied cell: index 81 tripped 8/8, indices 67 and
# 69 tripped 1/8, and three control indices tripped 0/8. A single attempt
# therefore misses real mode-A failures; the rate over several does not.
#
# Production makes 10 folds x 3 estimands x 2 arms = 60 such fits per run, so
# even a 1-in-8 per-fit rate on one bad cell is near-certain death - which is
# why a small nonzero sl_fail_rate here is a positive finding, not noise.
N_PROBE_CELLS <- 3L
SL_ATTEMPTS <- 3L

# ---- arguments --------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[[1]])
}
arg_ints <- function(name, default = NULL) {
  raw <- arg_value(name)
  if (is.null(raw)) return(default)
  as.integer(strsplit(raw, ",", fixed = TRUE)[[1]])
}

stages <- arg_ints("stage", c(1L, 2L))
n_controls <- as.integer(arg_value("controls", "100"))
out_dir <- arg_value("out", file.path(study$res_path, "diagnostics"))
log_dir <- arg_value("logs")
stage3_ids <- arg_ints("ids", c(67L, 68L))

failed_idx <- as.integer(readLines(study$failed_file))
failed_idx <- failed_idx[!is.na(failed_idx)]

# Keep the diagnostic single-threaded and deterministic. The production runs use
# plan(multisession, workers = 2); the failure does not depend on that, and
# sequential keeps the per-arm tracebacks in stage 3 readable.
plan(sequential)

cat(sprintf(
  paste0(
    "grf %s | pseudo %s | SuperLearner %s | %s\n",
    "failed indices: %d, from %s\n\n"
  ),
  as.character(packageVersion("grf")),
  as.character(packageVersion("pseudo")),
  as.character(packageVersion("SuperLearner")),
  R.version.string,
  length(failed_idx),
  study$failed_file
))

# ---- shared helpers ---------------------------------------------------------

#' The dataset one array index produced, exactly as surv_analysis.R built it
#'
#' Same setup_rng_stream(run) and same generate_surv_data() call. return_truth is
#' FALSE because nothing here scores against the truth and the per-individual
#' integrate() calls are the slowest part of the DGM.
index_data <- function(i, return_truth = FALSE) {
  param <- study$grid[i, ]
  setup_rng_stream(param$run)
  gen <- generate_surv_data(
    scenario = param$scenario,
    n = param$n,
    return_truth = return_truth,
    censoring = param$censoring
  )
  d <- gen$dataset
  list(
    param = param,
    X = as.matrix(d[, !(names(d) %in% c("Y", "D", "W"))]),
    Y = d$Y,
    D = d$D,
    W = d$W,
    # surv_analysis.R:44
    n_folds = ifelse(param$n < 300, 5, 10)
  )
}

#' The fold labels all_cate_surv_models() would build for n rows
#'
#' surv_models.R:30. Deterministic, so it needs no RNG state to match.
fold_labels <- function(n_obs, n_folds) sort(seq(n_obs) %% n_folds) + 1

#' Base R gained %||% in 4.4.0; the cluster runs 4.3.2, so define it either way
`%||%` <- function(x, y) if (is.null(x)) y else x

# =============================================================================
if (1L %in% stages) {
  cat("=== 1. What the failed indices are ===\n\n")

  failed_grid <- study$grid[failed_idx, ]

  by_scenario <- study$grid %>%
    mutate(failed = seq_len(n()) %in% failed_idx) %>%
    group_by(scenario) %>%
    summarise(cells = n(), failed = sum(failed), .groups = "drop") %>%
    left_join(
      survival_scenario_params[, c(
        "scenario", "description", "bW_1", "bW_2", "b3_1", "b3_2"
      )],
      by = "scenario"
    )

  cat("failures by scenario, against the DGM effects that define it:\n")
  print(as.data.frame(by_scenario), row.names = FALSE)
  cat(
    "\n  bW_1 is the treatment effect on the log scale of the EVENT OF INTEREST.\n",
    " The scenarios with bW_1 == 0 are the ones that never fail.\n\n"
  )

  cat("failures by censoring:\n")
  print(table(censoring = failed_grid$censoring))

  runs <- table(failed_grid$run)
  cat(sprintf(
    paste0(
      "\nfailures by run: %d distinct runs affected, %d..%d failures each.\n",
      "  10 = all five failing scenarios x both censoring settings. Clustering by\n",
      "  run is expected: setup_rng_stream(run) gives every scenario the same\n",
      "  draws of W, X1, X2, X3 and u, so a bad draw hits all of them at once.\n"
    ),
    length(runs), min(runs), max(runs)
  ))
  print(table(failures_in_that_run = as.integer(runs)))
  cat("\n")
}

# =============================================================================
# Cheap probe: everything that predicts the failure, with no model fitting
# beyond one SuperLearner fit on the single worst cell.
probe_index <- function(i) {
  d <- index_data(i)
  Y <- d$Y
  D_int <- as.integer(d$D)
  Dc <- as.integer(d$D %in% c(1, 2))
  W <- d$W
  n_obs <- length(Y)

  fold_indices <- fold_labels(n_obs, d$n_folds)
  fold_list <- unique(fold_indices)
  n_folds <- length(fold_list)

  # The whole-sample and leave-one-fold-out pseudo-values the guarded arms use,
  # so the README's ci.omit NA candidate is measured rather than assumed.
  ps_whole <- pseudo_all(Y, d$D, HORIZON)
  na_whole <- sum(is.na(unlist(ps_whole)))

  # Replicate the 3-way split of pseudo_sl_t_split (surv_models.R:806-821) and
  # look at the training-fold pseudo-values it hands to SuperLearner.
  cells <- list()
  na_split_train <- 0L
  for (k in seq_along(fold_list)) {
    km_fold <- fold_list[(k %% n_folds) + 1]
    in_train <- fold_indices != fold_list[k] & fold_indices != km_fold

    ps_yl <- pseudoyl(Y[in_train], D_int[in_train], HORIZON)
    ps_c <- pseudomean(Y[in_train], Dc[in_train], HORIZON)
    na_split_train <- na_split_train +
      sum(is.na(ps_yl$pseudo$cause1)) +
      sum(is.na(ps_yl$pseudo$cause2)) +
      sum(is.na(ps_c))

    W_train <- W[in_train]
    for (e in ESTIMANDS) {
      ps <- switch(
        e,
        RMTL1 = ps_yl$pseudo$cause1,
        RMTL2 = ps_yl$pseudo$cause2,
        RMSTc = ps_c
      )
      for (w in c(0L, 1L)) {
        y <- ps[W_train == w]
        cells[[length(cells) + 1]] <- data.frame(
          fold = fold_list[k],
          estimand = e,
          arm = w,
          n = length(y),
          sd = sd(y),
          n_unique = length(unique(round(y, 10))),
          # SuperLearner CV-folds a constant y whenever one value dominates
          max_tie_frac = max(table(round(y, 10))) / length(y)
        )
      }
    }
  }
  cells <- bind_rows(cells)

  # Reproduce the crash condition itself on the least varied cells: a NULL entry
  # in $fitLibrary is exactly what make_t_cate()'s predict(..., onlySL = FALSE)
  # then chokes on. See the SL_ATTEMPTS comment for why each cell is fitted more
  # than once. min_n_unique below is the continuous risk score; this is the
  # measurement of how often that score actually bites.
  cells <- cells[order(cells$n_unique), ]
  worst <- cells[1, ]

  null_fits <- 0L
  attempts <- 0L
  trips <- 0L
  predict_ok <- TRUE
  sl_note <- "ok"
  sl_stderr <- ""

  for (r in seq_len(min(N_PROBE_CELLS, nrow(cells)))) {
    cell <- cells[r, ]
    km_fold <- fold_list[(which(fold_list == cell$fold) %% n_folds) + 1]
    in_train <- fold_indices != cell$fold & fold_indices != km_fold
    ps_yl_c <- pseudoyl(Y[in_train], D_int[in_train], HORIZON)
    ps <- switch(
      cell$estimand,
      RMTL1 = ps_yl_c$pseudo$cause1,
      RMTL2 = ps_yl_c$pseudo$cause2,
      RMSTc = pseudomean(Y[in_train], Dc[in_train], HORIZON)
    )
    keep_w <- W[in_train] == cell$arm
    y_cell <- ps[keep_w]
    x_cell <- as.data.frame(d$X[in_train, , drop = FALSE][keep_w, , drop = FALSE])

    for (attempt in seq_len(SL_ATTEMPTS)) {
      attempts <- attempts + 1L

      # SuperLearner cats its learners' failures straight to stderr; keep them
      # out of the report but hang on to the first, it names the learner.
      fit <- NULL
      msg <- capture.output(
        fit <- tryCatch(
          suppressWarnings(SuperLearner(
            Y = y_cell, X = x_cell,
            SL.library = DEFAULT_SL_LIBRARY, cvControl = list(V = 5)
          )),
          error = function(e) conditionMessage(e)
        ),
        type = "message"
      )
      if (length(msg) && !nzchar(sl_stderr)) {
        sl_stderr <- trimws(paste(utils::head(msg, 2), collapse = " "))
      }

      if (is.character(fit)) {
        trips <- trips + 1L
        null_fits <- NA_integer_
        predict_ok <- FALSE
        sl_note <- paste0("SuperLearner() errored on ", cell$estimand,
                          "/W", cell$arm, ": ", fit)
        worst <- cell
        next
      }

      n_null <- sum(vapply(fit$fitLibrary, is.null, logical(1)))
      if (n_null == 0) next

      trips <- trips + 1L
      if (isTRUE(n_null > null_fits)) null_fits <- n_null
      worst <- cell

      pred <- tryCatch(
        suppressWarnings(predict(fit, newdata = x_cell[1:2, , drop = FALSE])),
        error = function(e) conditionMessage(e)
      )
      if (is.character(pred)) {
        predict_ok <- FALSE
        sl_note <- paste0("predict() errored on ", cell$estimand,
                          "/W", cell$arm, ": ", pred)
      }
    }
  }

  data.frame(
    index = i,
    scenario = d$param$scenario,
    censoring = d$param$censoring,
    run = d$param$run,
    # the driver: how much of the competing event survives to the horizon,
    # in each treatment arm
    ce_before_horizon_W0 = sum(d$D == 2 & Y <= HORIZON & W == 0),
    ce_before_horizon_W1 = sum(d$D == 2 & Y <= HORIZON & W == 1),
    at_risk_past_horizon = sum(Y > HORIZON),
    na_pseudo_whole = na_whole,
    na_pseudo_split_train = na_split_train,
    # the risk score: the fewest distinct pseudo-values any (fold, estimand,
    # treatment arm) cell hands to SuperLearner
    min_n_unique = cells$n_unique[1],
    min_estimand = cells$estimand[1],
    min_arm = cells$arm[1],
    max_tie_frac = round(max(cells$max_tie_frac), 4),
    # the cell that actually tripped, where one did
    worst_estimand = worst$estimand,
    worst_arm = worst$arm,
    worst_fold = worst$fold,
    sl_null_fits = null_fits,
    sl_attempts = attempts,
    sl_trips = trips,
    sl_fail_rate = round(trips / attempts, 3),
    sl_predict_ok = predict_ok,
    sl_note = sl_note,
    sl_stderr = sl_stderr,
    # Which of the two routes into pseudo_sl_t_split()'s SuperLearner() call
    # this index is set up for. Mode B is checked first: an NA aborts the fit
    # outright, so it wins wherever both are present.
    predicted_mode = if (na_split_train > 0) {
      "B_pseudo_NA"
    } else if (trips > 0) {
      "A_degenerate_library"
    } else {
      "unexplained"
    },
    stringsAsFactors = FALSE
  )
}

if (2L %in% stages) {
  cat("=== 2. Degeneracy probe over failed indices and controls ===\n\n")

  # --ids= narrows the failed side, for a quick smoke test of the probe itself.
  # Left unset, stage 2 sweeps every index in failed_ids.txt.
  probe_failed <- if (is.null(arg_value("ids"))) failed_idx else stage3_ids

  succeeded <- setdiff(seq_len(nrow(study$grid)), failed_idx)
  set.seed(1998) # controls are a fixed sample, so reruns are comparable
  controls <- sort(sample(succeeded, min(n_controls, length(succeeded))))

  todo <- c(probe_failed, controls)
  labels <- rep(c("failed", "control"), c(length(probe_failed), length(controls)))
  cat(sprintf(
    "probing %d failed + %d control indices (no model fitting except one\nSuperLearner fit per index, on its least varied cell)\n\n",
    length(probe_failed), length(controls)
  ))

  rows <- vector("list", length(todo))
  for (j in seq_along(todo)) {
    t0 <- Sys.time()
    rows[[j]] <- tryCatch(
      cbind(group = labels[j], probe_index(todo[j])),
      error = function(e) {
        data.frame(
          group = labels[j], index = todo[j], sl_note = paste0("PROBE ERROR: ", conditionMessage(e)),
          stringsAsFactors = FALSE
        )
      }
    )
    cat(sprintf(
      "  [%3d/%3d] index %4d (%-7s) %5.1fs  min_unique=%-3s sl_trips=%s/%s pseudo_NA=%-3s -> %s\n",
      j, length(todo), todo[j], labels[j],
      as.numeric(difftime(Sys.time(), t0, units = "secs")),
      rows[[j]]$min_n_unique %||% NA,
      rows[[j]]$sl_trips %||% NA,
      rows[[j]]$sl_attempts %||% NA,
      (rows[[j]]$na_pseudo_whole %||% NA) + (rows[[j]]$na_pseudo_split_train %||% NA),
      rows[[j]]$predicted_mode %||% "probe error"
    ))
    flush.console()
  }

  probe <- bind_rows(rows)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, "failed_probe.csv")
  write.csv(probe, out_file, row.names = FALSE)
  cat(sprintf("\nper-index probe written to %s\n\n", out_file))

  cat("does the SuperLearner library degenerate? (failed vs control)\n")
  print(probe %>%
    group_by(group) %>%
    summarise(
      n = n(),
      any_sl_trip = sum(sl_trips > 0, na.rm = TRUE),
      mean_sl_fail_rate = round(mean(sl_fail_rate, na.rm = TRUE), 3),
      predict_fails = sum(!sl_predict_ok, na.rm = TRUE),
      median_min_unique = median(min_n_unique, na.rm = TRUE),
      median_ce_W1 = median(ce_before_horizon_W1, na.rm = TRUE),
      median_ce_W0 = median(ce_before_horizon_W0, na.rm = TRUE),
      any_pseudo_NA = sum(na_pseudo_whole + na_pseudo_split_train > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>% as.data.frame(), row.names = FALSE)

  cat("\nwhich cell is the least varied one, among failed indices:\n")
  print(with(
    probe[probe$group == "failed", ],
    table(estimand = min_estimand, treatment_arm = min_arm)
  ))

  cat("\nrisk score (min_n_unique) by group:\n")
  print(with(probe, table(
    min_n_unique = cut(min_n_unique, c(0, 3, 5, 10, 20, Inf),
                       labels = c("<=3", "4-5", "6-10", "11-20", ">20")),
    group = group
  )))
  cat("\npredicted route into pseudo_sl_t_split()'s SuperLearner() call:\n")
  print(with(probe, table(predicted_mode = predicted_mode, group = group)))
  cat(
    "\n  A_degenerate_library  near-constant pseudo-values leave NULL entries in\n",
    "                       $fitLibrary, and predict(onlySL = FALSE) hits one\n",
    " B_pseudo_NA           unguarded pseudoyl/pseudomean returned NA, which\n",
    "                       SuperLearner() rejects outright\n",
    " unexplained           neither - stage 3 needs pointing at these\n"
  )

  # Say plainly how much of the failed list this probe accounts for. Do not read
  # silence as agreement: an unexplained failure is a third mode until shown
  # otherwise.
  unexplained <- probe[probe$group == "failed" &
                         probe$predicted_mode == "unexplained", ]
  false_alarms <- probe[probe$group == "control" &
                          probe$predicted_mode != "unexplained", ]
  cat(sprintf(
    paste0(
      "\n  %d of %d failed indices unexplained; %d of %d controls flagged anyway.\n",
      "  Run stage 3 on the unexplained ones before trusting the story above:\n",
      "    --stage=3 --ids=%s\n\n"
    ),
    nrow(unexplained), sum(probe$group == "failed"),
    nrow(false_alarms), sum(probe$group == "control"),
    paste(utils::head(unexplained$index, 6), collapse = ",")
  ))
}

# =============================================================================
#' Every arm of all_cate_surv_models(), each one wrapped so a failure is
#' recorded instead of aborting
#'
#' The point of not reusing all_cate_surv_models() directly is that it stops at
#' the first error, so a run that fails in two places only ever reports one. The
#' call order and arguments mirror surv_models.R:36-175.
debug_arms <- function(i) {
  d <- index_data(i)
  X <- d$X
  Y <- d$Y
  D <- d$D
  W <- d$W
  fold_indices <- fold_labels(length(Y), d$n_folds)
  fold_list <- unique(fold_indices)

  rows <- list()
  run_arm <- function(label, expr) {
    t0 <- Sys.time()
    tb <- NULL
    res <- tryCatch(
      withCallingHandlers(
        {
          value <- force(expr)
          list(status = "ok", message = sprintf("NA predictions: %d", sum(is.na(unlist(value)))))
        },
        # A calling handler sees the stack before it unwinds, which tryCatch's
        # handler does not. Keep first lines only and drop this harness's own
        # frames, or the useful part is buried under tryCatch machinery.
        error = function(e) {
          calls <- vapply(sys.calls(), function(cl) deparse(cl)[1], character(1))
          calls <- calls[!grepl(
            paste0(
              "^(tryCatch|tryCatchList|tryCatchOne|doTryCatch|withCallingHandlers",
              "|force|run_arm|debug_arms|lapply|bind_rows|FUN|function",
              "|\\.handleSimpleError|h\\(|stop\\()"
            ),
            calls
          )]
          tb <<- paste(utils::tail(calls, 5), collapse = " -> ")
        }
      ),
      error = function(e) list(status = "ERROR", message = conditionMessage(e))
    )
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    msg <- gsub("[\r\n]+", " ", res$message)
    cat(sprintf("  %-32s %7.1fs  %-5s %s\n", label, secs, res$status, msg))
    if (res$status == "ERROR" && !is.null(tb)) {
      cat(sprintf("      calls: %s\n", tb))
    }
    flush.console()
    rows[[length(rows) + 1]] <<- data.frame(
      index = i, arm = label, seconds = round(secs, 1),
      status = res$status, message = msg, traceback = tb %||% "",
      stringsAsFactors = FALSE
    )
    invisible(NULL)
  }

  cat(sprintf(
    "== index %d: scenario %d, n %d, censoring %s, run %d\n",
    i, d$param$scenario, d$param$n, d$param$censoring, d$param$run
  ))

  for (e in c(1, 2, "composite")) {
    run_arm(paste0("ipw ", e), cf_ipw(X, Y, D, W, HORIZON, event = e))
  }
  for (e in c(1, 2, "composite")) {
    run_arm(paste0("csf_cs ", e), csf_cs(X, Y, D, W, HORIZON, event = e))
  }
  for (e in c(1, 2)) {
    run_arm(paste0("csf_sh ", e), csf_sh(X, Y, D, W, HORIZON, event = e))
  }

  ps_whole <- pseudo_all(Y, D, HORIZON)
  ps_cv <- pseudo_crossfit(Y, D, HORIZON, fold_indices, fold_list, ps_whole)
  cat(sprintf(
    "  %-32s          pseudo NAs: whole %d, cvps n_na_fallback %s\n",
    "(pseudo-values)", sum(is.na(unlist(ps_whole))),
    paste(ps_cv$n_na_fallback, collapse = "/")
  ))

  pw <- function(e) ps_whole[[paste0("ps_", e)]]
  pc <- function(e) ps_cv[[paste0("ps_", e)]]

  for (e in ESTIMANDS) {
    run_arm(paste0("pseudo_cf_whole_oob ", e), pseudo_cf_whole_oob(X, pw(e), W))
    run_arm(paste0("pseudo_cf_whole_scf ", e),
            pseudo_cf_scf(X, pw(e), W, fold_indices, fold_list))
    run_arm(paste0("pseudo_cf_cvps_scf ", e),
            pseudo_cf_scf(X, pc(e), W, fold_indices, fold_list))
  }

  for (e in ESTIMANDS) {
    nuis_oob <- NULL
    run_arm(paste0("nuis_rf_whole_oob ", e), {
      nuis_oob <- nuisance_pseudo_rf_oob(X, pw(e), W)
      nuis_oob$po
    })
    if (!is.null(nuis_oob)) {
      run_arm(paste0("dr_whole_oob ", e), stage2_whole_rf(X, nuis_oob$po)$tau)
    }
    nuis_scf <- NULL
    run_arm(paste0("nuis_rf_cvps_scf ", e), {
      nuis_scf <- nuisance_pseudo_rf_scf(X, pc(e), W, pw(e), fold_indices, fold_list)
      nuis_scf$po
    })
    if (!is.null(nuis_scf)) {
      run_arm(paste0("dr_cvps_scf ", e),
              stage_2_rf_scf(X, nuis_scf$po, fold_indices, fold_list))
    }
  }

  for (e in ESTIMANDS) {
    run_arm(paste0("sl_t_whole ", e),
            pseudo_sl_t_standard(X, pw(e), W, fold_indices, fold_list, DEFAULT_SL_LIBRARY))
    run_arm(paste0("sl_t_cvps ", e),
            pseudo_sl_t_standard(X, pc(e), W, fold_indices, fold_list, DEFAULT_SL_LIBRARY))
    nuis_sl <- NULL
    run_arm(paste0("nuis_sl_whole ", e), {
      nuis_sl <- nuisance_pseudo_sl(X, pw(e), W, pw(e), fold_indices, fold_list,
                                    DEFAULT_SL_LIBRARY)
      nuis_sl$po
    })
    if (!is.null(nuis_sl)) {
      run_arm(paste0("sl_dr_whole ", e),
              pseudo_dr_sl(X, nuis_sl$po, fold_indices, fold_list, DEFAULT_SL_LIBRARY))
    }
  }

  # The arm all 225 failures came from. Kept last so everything else is already
  # reported when it goes. pseudo_whole is passed because all_cate_surv_models()
  # passes it - without it this would exercise the standalone drop-the-NA-row
  # path rather than the production one.
  run_arm("sl_t_split (all estimands)",
          pseudo_sl_t_split(X, Y, D, W, HORIZON, fold_indices, fold_list,
                            DEFAULT_SL_LIBRARY, pseudo_whole = ps_whole))

  cat("\n")
  bind_rows(rows)
}

if (3L %in% stages) {
  cat("=== 3. Per-arm rerun, continuing past failures ===\n\n")
  cat(sprintf("indices: %s\n\n", paste(stage3_ids, collapse = ", ")))

  arm_rows <- bind_rows(lapply(stage3_ids, debug_arms))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, "failed_arms.csv")
  write.csv(arm_rows, out_file, row.names = FALSE)
  cat(sprintf("per-arm results written to %s\n\n", out_file))

  failures <- arm_rows[arm_rows$status == "ERROR", ]
  if (nrow(failures) == 0) {
    cat("  No arm errored. If these indices really did fail on the cluster,\n")
    cat("  the cause is environmental (R 4.3.2 there vs 4.5.3 here) - go to stage 4.\n\n")
  } else {
    cat("arms that errored:\n")
    print(as.data.frame(failures[, c("index", "arm", "message")]), row.names = FALSE)
    cat("\n")
  }
}

# =============================================================================
if (4L %in% stages) {
  cat("=== 4. PBS .o log triage ===\n\n")

  if (is.null(log_dir) || !dir.exists(log_dir)) {
    stop("stage 4 needs --logs=<directory of PBS .o files>")
  }

  # all_cate_surv_models()'s progress message()s go to stderr, so a .o file only
  # carries surv_analysis.R's three print() calls. That is still enough to say
  # whether the job died in the models, at saveRDS, or never ran.
  files <- list.files(log_dir, pattern = "\\.o[0-9]", full.names = TRUE)
  cat(sprintf("%d log files in %s\n\n", length(files), log_dir))

  # PBS array logs are <jobname>.o<jobid>.<array index> (or ...[index]). A
  # non-array log has no trailing index, and the regex would then carve one out
  # of the job id, so anything outside the grid's row range is dropped rather
  # than quietly counted against the wrong index.
  classify <- function(f) {
    txt <- tryCatch(readLines(f, warn = FALSE), error = function(e) character())
    idx <- suppressWarnings(as.integer(sub(".*\\.o[0-9]+[.\\[]?([0-9]+).*", "\\1", basename(f))))
    if (!is.na(idx) && (idx < 1 || idx > nrow(study$grid))) idx <- NA_integer_
    # surv_analysis.R's options(error=) handler writes the reason to stdout, so
    # for anything run after that was added the log names its own cause and
    # there is no need to reproduce it locally at all.
    reported <- grep("^error:", txt, value = TRUE)

    verdict <- if (!length(txt)) {
      "empty"
    } else if (any(grepl("=== SIMULATION FAILED ===", txt, fixed = TRUE))) {
      "reported_error"
    } else if (any(grepl("Simulation completed!", txt, fixed = TRUE))) {
      "completed"
    } else if (any(grepl("Time difference of", txt, fixed = TRUE))) {
      "died_at_save"
    } else if (any(grepl("scenario", txt))) {
      "died_in_models"
    } else {
      "unrecognised"
    }
    data.frame(file = basename(f), array_index = idx, verdict = verdict,
               error = if (length(reported)) reported[[1]] else "",
               n_lines = length(txt), stringsAsFactors = FALSE)
  }

  logs <- bind_rows(lapply(files, classify))

  # surv_rerun.sh's PBS_ARRAY_INDEX is a LINE NUMBER in failed_ids.txt, not a
  # grid row - it seds the grid index out of the file. Pass --rerun for logs
  # from logs_rerun/ so the two numbering schemes do not get mixed up.
  if ("--rerun" %in% args) {
    cat("--rerun: treating the array index in each filename as a line of failed_ids.txt\n\n")
    logs$line_in_failed_ids <- logs$array_index
    logs$array_index <- failed_idx[logs$line_in_failed_ids]
  }
  logs$in_failed_ids <- logs$array_index %in% failed_idx

  cat("verdict by whether the index is in failed_ids.txt:\n")
  print(table(verdict = logs$verdict, in_failed_ids = logs$in_failed_ids))
  cat(
    "\n  reported_error  surv_analysis.R's error handler named the cause itself -\n",
    "                 see the `error` column, no local reproduction needed\n",
    " died_in_models  only the param row printed, and no handler output: either a\n",
    "                 log from before the handler existed, or a scheduler kill\n",
    " died_at_save    the t1-t0 timing printed but not the completion message,\n",
    "                 so saveRDS failed: disk or quota, not the model code\n",
    " completed       the RDS was written; check_failed is looking in the wrong place\n",
    " empty           the job never started, or was killed before any output\n"
  )

  if (any(nzchar(logs$error))) {
    cat("\nerrors the jobs reported for themselves:\n")
    print(as.data.frame(table(error = logs$error[nzchar(logs$error)])),
          row.names = FALSE)
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, "log_triage.csv")
  write.csv(logs, out_file, row.names = FALSE)
  cat(sprintf("\nper-log verdicts written to %s\n\n", out_file))
}
