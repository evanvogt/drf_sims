##########
# Title: why the array indices in jobscripts/failed_ids.txt fail
##########
# The HPC is only keeping the PBS `.o` files, so the R error that killed these
# runs is invisible on the cluster. This script reproduces it locally. It
# diagnoses only - it changes nothing in surv_models.R.
#
# STATUS: two rounds, two different answers. Read both.
#
# ROUND 1 (225 failures of 1400) is written up below and is FIXED - see
# pseudo_sl_t_split()'s header comment in surv_models.R.
#
# ROUND 2 (19 failures) is what stages 1-3 now measure. Those 19 are all a
# SUBSET of the original 225, and they came back unchanged from a rerun at 4h
# and 20gb, so they are deterministic in the array index and not a resource
# problem. They did not fail for the reasons round 1 failed, and the mode labels
# this script used in round 1 would have hidden that - so they have been
# rewritten. Do not restore the old ones.
#
#   MODE C - whole-sample pseudo-value NA (15 of the 19, censoring = TRUE)
#
#   The round-1 fix guarded pseudo_sl_t_split()'s TRAINING-fold pseudo-values.
#   pseudo_all() itself is still unguarded - it is the one pseudo-value producer
#   in surv_models.R with no is.na() fallback - and pseudo_crossfit() fills its
#   own NAs FROM pseudo_all()'s vector, so when that vector is NA the fallback
#   substitutes an NA and the cvps arms inherit it too.
#
#   It goes NA when nobody is observed past the horizon (at_risk_past_horizon
#   == 0): pseudoyl()'s leave-one-out risk set then empties at or below tmax, so
#   the pseudo::ci.omit bug in README "Known issues" reaches the estimand rather
#   than the unused tail. Confirmed on index 295 - max(Y) = 27.83 against
#   horizon 28, one NA on each of cause 1 and cause 2.
#
#   The first consumer is pseudo_cf_whole_oob(), which is the FIRST pseudo-value
#   arm in all_cate_surv_models() and long before any SuperLearner one, and grf
#   refuses an NA outcome:
#
#     The vector of observations (W, Y, Z or D) contains at least one NA.
#
#   So blaming a SuperLearner arm for these would be blaming an arm the run
#   never reaches. cf_whole_na_error in stage 2 is the receipt.
#
#   MODE A' - pretest's own fallback fails (4 of the 19, censoring = FALSE)
#
#   Round 1's mode A condition, but NOT round 1's mode A failure, and this is the
#   distinction that matters: the round-1 fix added pretest_superlearner() and
#   onlySL = TRUE, and neither of them covers this.
#
#   These four are more extreme than anything round 1 saw. min_n_unique is 1 - a
#   LITERALLY constant treated-arm RMTL2 cell - against 3 on index 67. On a
#   constant y, pretest_superlearner() drops every real candidate and takes its
#   own "every candidate failed/warned" branch (R/cate_models.R:415-425), which
#   returns SL.mean. The live SuperLearner() call then drops SL.mean too and
#   errors before there is any fit at all:
#
#     All algorithms dropped from library
#
#   Confirmed on index 923: sl_pretest_kept = "RMTL2/W1 kept: SL.mean", followed
#   by that error on 3 of 9 attempts. onlySL = TRUE is irrelevant here - the
#   failure is at fit time, not predict time - so this is a live route in EVERY
#   SuperLearner arm including the "fixed" pseudo_sl_t_split().
#
#   Stage 2 separates the two by whether the fit or the predict died
#   (A_sl_fit_error vs A_null_fit_predict) and records both prediction paths
#   (sl_predict_ok, sl_predict_onlysl_ok), so which fix is needed is measured
#   rather than argued.
#
# ROUND 1 ROOT CAUSE (of the 225 failures of 1400, as found): every failure reproduced
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
# Stage 4 had NO INPUT for round 2: both the .e and the .o files for that rerun
# were deleted before anyone looked at them, so surv_analysis.R's error handler -
# which writes the message and call stack to stdout precisely so a lost .e does
# not lose the cause - had nothing left to be read from. That is why round 2 is
# entirely local reproduction. Keep logs_rerun/ next time and most of stages 2
# and 3 becomes unnecessary.
#
# Costs, measured on R 4.5.3 locally with plan(sequential):
#   stage 1  instant
#   stage 2  ~12 s an index, so ~65 min for the 225 failed plus 100 controls,
#            or ~25 min for round 2's 19 plus 100 controls (most of it the
#            N_PROBE_CELLS x SL_ATTEMPTS SuperLearner fits)
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
  #
  # na_whole is the one that is still FATAL. pseudo_all() has no NA guard - it is
  # the only pseudo-value producer in surv_models.R without one - and
  # pseudo_crossfit()'s is.na() fallback substitutes from this same vector, so an
  # NA here reaches every pseudo-value arm rather than being filled. The split
  # arm's own NAs (na_split_train, below) are guarded as of the SL fix and are
  # now informational only. Do not let the two share a mode label.
  ps_whole <- pseudo_all(Y, d$D, HORIZON)
  na_whole <- sum(is.na(unlist(ps_whole)))

  # The assertion the whole-sample-NA route rests on: all_cate_surv_models()
  # reaches pseudo_cf_whole_oob() before any SuperLearner arm, and grf refuses an
  # NA outcome ("The vector of observations (W, Y, Z or D) contains at least one
  # NA."). Measured on the real vector rather than asserted from the docs, and
  # only when there is an NA to feed it.
  cf_whole_na_error <- ""
  if (na_whole > 0) {
    bad <- names(ps_whole)[vapply(ps_whole, anyNA, logical(1))]
    err <- tryCatch(
      {
        pseudo_cf_whole_oob(d$X, ps_whole[[bad[[1]]]], W)
        "NO ERROR - grf accepted the NA"
      },
      error = function(e) conditionMessage(e)
    )
    cf_whole_na_error <- paste0(paste(bad, collapse = "+"), ": ", err)
  }

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
  # in $fitLibrary is exactly what a predict(..., onlySL = FALSE) then chokes on.
  # See the SL_ATTEMPTS comment for why each cell is fitted more than once.
  # min_n_unique below is the continuous risk score; this is the measurement of
  # how often that score actually bites.
  #
  # The library is PRETESTED here, because production now pretests everywhere -
  # measuring the raw DEFAULT_SL_LIBRARY would re-measure the bug that was
  # already fixed. What is NOT fixed everywhere is the prediction call:
  # onlySL = TRUE was added to pseudo_sl_t_split() alone (surv_models.R:902),
  # while pseudo_sl_t_standard() (surv_models.R:775-776), nuisance_pseudo_sl()
  # and R/cate_models.R::stage_2_sl still take the onlySL = FALSE default. So a
  # candidate that survives pretest's 2-fold CV and then fails in the live 5-fold
  # one still leaves a NULL fit for those arms to trip on. Both prediction paths
  # are therefore recorded: predict_ok is the live one in the unfixed arms,
  # predict_onlysl_ok is what those arms would do if they were fixed like the
  # split arm. predict_ok FALSE with predict_onlysl_ok TRUE names the fix.
  cells <- cells[order(cells$n_unique), ]
  worst <- cells[1, ]

  null_fits <- 0L
  attempts <- 0L
  trips <- 0L
  predict_ok <- TRUE
  predict_onlysl_ok <- TRUE
  sl_note <- "ok"
  sl_stderr <- ""
  pretest_note <- ""

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

    # pretest_superlearner() cats its own removals to stdout; keep them out of
    # the report, but record what it kept - an empty-ish library is itself a
    # finding, and "SL.mean" means every candidate died on this cell.
    lib <- NULL
    capture.output(
      lib <- suppressWarnings(
        pretest_superlearner(y_cell, x_cell, DEFAULT_SL_LIBRARY, gaussian())
      )
    )
    if (!nzchar(pretest_note) && !identical(sort(lib), sort(DEFAULT_SL_LIBRARY))) {
      pretest_note <- paste0(
        cell$estimand, "/W", cell$arm, " kept: ", paste(lib, collapse = "+")
      )
    }

    for (attempt in seq_len(SL_ATTEMPTS)) {
      attempts <- attempts + 1L

      # SuperLearner cats its learners' failures straight to stderr; keep them
      # out of the report but hang on to the first, it names the learner.
      fit <- NULL
      msg <- capture.output(
        fit <- tryCatch(
          suppressWarnings(SuperLearner(
            Y = y_cell, X = x_cell,
            SL.library = lib, cvControl = list(V = 5)
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
        predict_onlysl_ok <- FALSE
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

      # The unfixed arms' call: default onlySL = FALSE predicts from every
      # library member, NULL fits included.
      pred <- tryCatch(
        suppressWarnings(predict(fit, newdata = x_cell[1:2, , drop = FALSE])),
        error = function(e) conditionMessage(e)
      )
      if (is.character(pred)) {
        predict_ok <- FALSE
        sl_note <- paste0("predict() errored on ", cell$estimand,
                          "/W", cell$arm, ": ", pred)
      }

      # The split arm's call, for contrast.
      pred_only <- tryCatch(
        suppressWarnings(predict(fit, newdata = x_cell[1:2, , drop = FALSE],
                                 onlySL = TRUE)),
        error = function(e) conditionMessage(e)
      )
      if (is.character(pred_only)) {
        predict_onlysl_ok <- FALSE
        sl_note <- paste0(sl_note, " | onlySL = TRUE errored too: ", pred_only)
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
    # the driver of the whole-sample NA route. pseudoyl()'s leave-one-out risk
    # set empties at the tail, so the max-time individual gets NaN - and that
    # only reaches the estimand when the max observed time is at or below the
    # horizon, i.e. exactly when at_risk_past_horizon is 0.
    at_risk_past_horizon = sum(Y > HORIZON),
    max_Y = round(max(Y), 2),
    na_pseudo_whole = na_whole,
    cf_whole_na_error = cf_whole_na_error,
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
    sl_predict_onlysl_ok = predict_onlysl_ok,
    sl_pretest_kept = pretest_note,
    sl_note = sl_note,
    sl_stderr = sl_stderr,
    # Which route this index is set up for, AGAINST THE CODE AS IT STANDS.
    #
    # This is not the taxonomy the 225-failure sweep used, and the difference
    # matters. Back then both routes ran through pseudo_sl_t_split() and mode B
    # was keyed on na_split_train. The SL fix guarded that arm, so na_split_train
    # alone is no longer a cause of anything - keying on it would relabel the
    # whole-sample NA route as "the bug we already fixed" and hide it.
    #
    # Checked in the order all_cate_surv_models() would hit them, so the label
    # names the arm that actually aborts the run first.
    predicted_mode = if (na_whole > 0) {
      # pseudo_all() is unguarded, and pseudo_crossfit() fills from it, so this
      # reaches pseudo_cf_whole_oob() - the FIRST pseudo-value arm, long before
      # any SuperLearner one. cf_whole_na_error is the receipt.
      "C_whole_pseudo_NA"
    } else if (trips > 0 && !predict_ok && !predict_onlysl_ok) {
      # The live SuperLearner() call itself died, so there was never a fit to
      # predict from. onlySL = TRUE cannot help. This is where pretest's own
      # SL.mean fallback lands on a constant cell: pretest drops every real
      # candidate, returns SL.mean, and SuperLearner then drops that too
      # ("All algorithms dropped from library"). Distinct from round 1's mode A,
      # and NOT covered by the round-1 fix - see sl_pretest_kept.
      "A_sl_fit_error"
    } else if (trips > 0 && !predict_ok) {
      # Round 1's mode A proper: the fit succeeded but left a NULL entry in
      # $fitLibrary, and the arms that still predict with onlySL = FALSE hit it.
      # predict_onlysl_ok is TRUE here by construction, so adding onlySL = TRUE
      # to those call sites is the whole fix.
      "A_null_fit_predict"
    } else if (trips > 0) {
      # NULL fits appeared but neither prediction path errored: a risk, not a
      # cause. Do not count it as an explanation.
      "A_degenerate_survivable"
    } else if (na_split_train > 0) {
      # The old mode B. Guarded since the SL fix - recorded so its absence as a
      # cause is visible, not so it can be blamed again.
      "B_split_pseudo_NA_guarded"
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
      "  [%3d/%3d] index %4d (%-7s) %5.1fs  min_unique=%-3s sl_trips=%s/%s NA_whole=%-3s NA_split=%-3s -> %s\n",
      j, length(todo), todo[j], labels[j],
      as.numeric(difftime(Sys.time(), t0, units = "secs")),
      rows[[j]]$min_n_unique %||% NA,
      rows[[j]]$sl_trips %||% NA,
      rows[[j]]$sl_attempts %||% NA,
      rows[[j]]$na_pseudo_whole %||% NA,
      rows[[j]]$na_pseudo_split_train %||% NA,
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
      # how many of those predict_fails onlySL = TRUE would have survived: the
      # size of the prize for fixing the three call sites that still omit it
      onlysl_would_fix = sum(!sl_predict_ok & sl_predict_onlysl_ok, na.rm = TRUE),
      median_min_unique = median(min_n_unique, na.rm = TRUE),
      median_ce_W1 = median(ce_before_horizon_W1, na.rm = TRUE),
      median_ce_W0 = median(ce_before_horizon_W0, na.rm = TRUE),
      .groups = "drop"
    ) %>% as.data.frame(), row.names = FALSE)

  cat("\nis the WHOLE-SAMPLE pseudo-value NA? (the unguarded route)\n")
  print(probe %>%
    group_by(group) %>%
    summarise(
      n = n(),
      na_whole = sum(na_pseudo_whole > 0, na.rm = TRUE),
      na_split_only = sum(na_pseudo_whole == 0 & na_pseudo_split_train > 0,
                          na.rm = TRUE),
      nobody_past_horizon = sum(at_risk_past_horizon == 0, na.rm = TRUE),
      median_max_Y = round(median(max_Y, na.rm = TRUE), 1),
      .groups = "drop"
    ) %>% as.data.frame(), row.names = FALSE)

  cf_err <- unique(probe$cf_whole_na_error[nzchar(probe$cf_whole_na_error)])
  if (length(cf_err)) {
    cat("\nwhat pseudo_cf_whole_oob() does with that NA:\n")
    for (e in cf_err) cat("  ", e, "\n", sep = "")
  }

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
  cat("\npredicted route, against the code as it stands:\n")
  print(with(probe, table(predicted_mode = predicted_mode, group = group)))
  cat(
    "\n  C_whole_pseudo_NA         pseudo_all() returned NA and nothing guards it.\n",
    "                            Kills pseudo_cf_whole_oob(), the first pseudo-\n",
    "                            value arm. Driver: at_risk_past_horizon == 0\n",
    " A_sl_fit_error            the live SuperLearner() call died outright. On a\n",
    "                            constant cell pretest keeps only SL.mean and that\n",
    "                            is dropped too. onlySL = TRUE cannot help\n",
    " A_null_fit_predict        round 1's mode A: fit succeeded but left a NULL in\n",
    "                            $fitLibrary, and an arm predicting with\n",
    "                            onlySL = FALSE hit it. onlySL = TRUE is the fix\n",
    " A_degenerate_survivable   NULL fits appeared but no prediction path errored:\n",
    "                            a risk, not a cause\n",
    " B_split_pseudo_NA_guarded the old mode B. Guarded since the SL fix, so it is\n",
    "                            reported and NOT counted as an explanation\n",
    " unexplained               none of the above - stage 3 needs pointing at these\n"
  )

  # Say plainly how much of the failed list this probe accounts for. Do not read
  # silence as agreement: an unexplained failure is a further mode until shown
  # otherwise. Only the two FATAL modes count as explanations - a guarded NA or a
  # survivable NULL fit is a description of the data, not a cause of death.
  FATAL_MODES <- c("C_whole_pseudo_NA", "A_sl_fit_error", "A_null_fit_predict")
  unexplained <- probe[probe$group == "failed" &
                         !(probe$predicted_mode %in% FATAL_MODES), ]
  false_alarms <- probe[probe$group == "control" &
                          probe$predicted_mode %in% FATAL_MODES, ]
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
    # whole_scf was missing from this harness even though all_cate_surv_models()
    # runs it (surv_models.R:88-98). It is the arm that isolates the pseudo-value
    # factor from the fitting factor, and it consumes the whole-sample vector, so
    # leaving it out would let a whole-sample NA look arm-specific.
    nuis_whole_scf <- NULL
    run_arm(paste0("nuis_rf_whole_scf ", e), {
      nuis_whole_scf <- nuisance_pseudo_rf_scf(X, pw(e), W, pw(e), fold_indices,
                                               fold_list)
      nuis_whole_scf$po
    })
    if (!is.null(nuis_whole_scf)) {
      run_arm(paste0("dr_whole_scf ", e),
              stage_2_rf_scf(X, nuis_whole_scf$po, fold_indices, fold_list))
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
    # sl_dr_cvps was missing too (surv_models.R:158-168). Same reason: it is live
    # in production and it consumes both pseudo-value vectors.
    nuis_sl_cv <- NULL
    run_arm(paste0("nuis_sl_cvps ", e), {
      nuis_sl_cv <- nuisance_pseudo_sl(X, pc(e), W, pw(e), fold_indices, fold_list,
                                       DEFAULT_SL_LIBRARY)
      nuis_sl_cv$po
    })
    if (!is.null(nuis_sl_cv)) {
      run_arm(paste0("sl_dr_cvps ", e),
              pseudo_dr_sl(X, nuis_sl_cv$po, fold_indices, fold_list,
                           DEFAULT_SL_LIBRARY))
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
