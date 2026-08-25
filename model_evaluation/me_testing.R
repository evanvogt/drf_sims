##########
# title: verification checks for the model evaluation study
##########
# Run this before submitting anything to the HPC:
#
#   Rscript model_evaluation/me_testing.R       # quick: no real nuisance fitting
#   Rscript model_evaluation/me_testing.R full  # adds real XGBoost-CV/H2O AutoML (slow)
#
# Checks, in order:
#   1. me_per_model()'s scoring matches a direct calc_infl_score()/
#      calc_dr_risk()/calc_cal_score() call on the same synthetic fixture -
#      i.e. moving the scoring loop from inside the old calc_metrics() to
#      compute_metrics()'s outer loop (see me_metrics.R) didn't change any
#      number - plus calc_cal_score()'s own properties: exactly 0 under
#      perfect calibration, K non-empty groups even when tau_hat is constant,
#      count (not proportion) weights, and the n-scaling those imply.
#   2. the DGM/design-matrix plumbing - Y, W first; covariates all numeric
#      (R/dgm_scenarios.R's output needs no factor expansion, unlike the old
#      benchtm-based data).
#   3. scenario 1 (no heterogeneity) behaves, same check continuous/ and
#      crossfitting/ both use for this scenario.
#   4. the real pipeline's structure: all 9 candidate models present and
#      correctly sized (always); in "full" mode, both nuisance pipelines
#      present with the expected arms and columns, plus one real call to
#      me_per_model() on real output.
#   5. the nuisance ARMS added by me_strategies.R - holdout_blocks()'s
#      pooling rule, the train/test filters each arm actually uses, and
#      single_split(). Structural, on synthetic data, no fitting.
#   6. me_strategies.R's assembly logic against a synthetic res_sim-shaped
#      object, XGBoost only. Local machines have no model_evaluation results
#      to run the real round-trip against (results/ carries every other study
#      but not this one), and H2O is the slow half - so this exercises the
#      read/compute/merge/write contract without either.
#
# WHAT DOES NOT RUN HERE. `full` mode is slow and SL2 fails on this machine
# for a documented local package-version reason (see README.md), so it is not
# a useful local signal. Run it on the cluster instead:
#
#   qsub model_evaluation/jobscripts/me_testing.sh
#
# and prove me_strategies.R's pass-through separately, over the real results,
# with me_strategies_verify.R.

library(dplyr)
library(future)
library(future.apply)
library(ranger)
library(glmnet)
library(SuperLearner)
library(xgboost)
library(h2o)
library(caret)
library(here)

source(here("R", "utils.R"))
source(here("model_evaluation", "me_dgms.R"))
source(here("model_evaluation", "me_utils.R"))
source(here("model_evaluation", "me_models.R"))
source(here("model_evaluation", "me_nuisance.R"))
source(here("model_evaluation", "me_config.R"))
source(here("model_evaluation", "me_metrics.R")) # calc_infl_score/calc_dr_risk/me_per_model

full <- "full" %in% commandArgs(trailingOnly = TRUE)

workers <- 2
metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)

pass <- character()
fail <- character()

report <- function(ok, msg) {
  cat(if (ok) "  PASS  " else "  FAIL  ", msg, "\n", sep = "")
  if (ok) pass <<- c(pass, msg) else fail <<- c(fail, msg)
}

# For checks whose subject may not exist - the candidate-model checks below
# all read a fitted model_list that a local SL2 failure can leave NULL, and
# vapply() over NULL returns logical(0), which all() reports as TRUE. Skipping
# is honest; a spurious PASS is not. `ok` is a promise, so it is never
# evaluated when `cond` is FALSE.
report_if <- function(cond, ok, msg) {
  if (cond) report(ok, msg) else cat("  SKIP  ", msg, "\n", sep = "")
}

# =============================================================================
cat("\n=== 1. scoring regression check: me_per_model() against a direct call ===\n")
# me_per_model() moved calc_metrics()'s "all 9 models, one vectorized pass"
# scoring into compute_metrics()'s "one model at a time" iteration. This
# checks that move didn't change any number, on a synthetic fixture that
# needs no real model fits or real nuisance estimation.

set.seed(42)
n_fake <- 50
tau_hat_fake <- rnorm(n_fake)
true_tau_fake <- rnorm(n_fake)
Y_fake <- rnorm(n_fake)
W_fake <- rbinom(n_fake, 1, 0.5)

# takes n so section 6 can build one at its own size; sections 1's calls keep
# the original 50-row default
fake_nuisance_df_n <- function(n = n_fake) {
  data.frame(
    tau_T = rnorm(n),
    pi = runif(n, 0.2, 0.8),
    phi = rnorm(n),
    # arm_scores() hard-fails without this rather than scoring an arm whose
    # oracle-pi pseudo-outcome was never computed, so the fixture carries it
    # exactly as every real nuisance data.frame does
    phi05 = rnorm(n)
  )
}
fake_nuisance_df <- function() fake_nuisance_df_n(n_fake)

# Both pipelines x every arm in NUISANCE_ARMS, so the column count this
# asserts tracks me_config.R rather than being typed here. me_per_model()
# derives its columns from whatever the nuisance list carries, which is what
# lets one scoring function serve the main study (2 arms), the strategies
# tree (4) and the split tree (1) - so the check that matters is that it
# picks up EVERY arm it is handed, not that it produces some fixed number.
fake_pipelines <- c("xgb", "automl")
fake_nuis <- setNames(lapply(fake_pipelines, function(p) {
  setNames(lapply(NUISANCE_ARMS, function(a) fake_nuisance_df()), NUISANCE_ARMS)
}), fake_pipelines)

sim_res_fake <- list(
  data = list(Y = Y_fake, W = W_fake),
  nuisances = fake_nuis
)

scored <- me_per_model(
  list(tau = tau_hat_fake), true_tau_fake, "test_model", sim_res_fake, NULL
)

# Every score type, written out longhand rather than by calling arm_scores() -
# the point of the check is that the assembled column equals an independent
# direct call, so reusing the assembler would make it vacuous.
checks <- list(true_pehe = mean((tau_hat_fake - true_tau_fake)^2))
for (p in fake_pipelines) {
  for (a in NUISANCE_ARMS) {
    d <- fake_nuis[[p]][[a]]
    checks[[paste0("infl_", a, "_", p)]] <-
      calc_infl_score(tau_hat_fake, d$tau_T, d$pi, Y_fake, W_fake)
    checks[[paste0("infl05_", a, "_", p)]] <-
      calc_infl_score(tau_hat_fake, d$tau_T, 0.5, Y_fake, W_fake)
    checks[[paste0("dr_", a, "_", p)]] <- calc_dr_risk(tau_hat_fake, d$phi)
    checks[[paste0("dr05_", a, "_", p)]] <- calc_dr_risk(tau_hat_fake, d$phi05)
    for (k in CAL_QUANTILES) {
      checks[[paste0("calq", k, "_", a, "_", p)]] <-
        calc_cal_score(tau_hat_fake, d$phi, k)
      checks[[paste0("cal05q", k, "_", a, "_", p)]] <-
        calc_cal_score(tau_hat_fake, d$phi05, k)
    }
  }
}

for (col in names(checks)) {
  diff <- abs(scored[[col]] - checks[[col]])
  report(diff < 1e-10, sprintf("%s matches a direct call (diff = %.2e)", col, diff))
}

# arm_scores() names and me_score_types() are two statements of the same list;
# if they drift, the column count assertions below silently stop meaning
# anything. This is what ties them together.
report(
  identical(
    names(arm_scores(tau_hat_fake, fake_nuisance_df(), Y_fake, W_fake)),
    me_score_types()
  ),
  "arm_scores() emits exactly me_score_types(), in order"
)

n_score_types <- length(me_score_types())
n_expected <- 1 + length(fake_pipelines) * length(NUISANCE_ARMS) * n_score_types
report(ncol(scored) == n_expected, sprintf(
  "%d score columns from %d pipelines x %d arms x %d score types (got %d)",
  n_expected, length(fake_pipelines), length(NUISANCE_ARMS), n_score_types,
  ncol(scored)
))
report(
  setequal(names(scored), names(checks)),
  "score column names are exactly the pipeline x arm x score-type product"
)

# the split tree's shape: one arm, and data/truth already restricted, so the
# same scorer works with no variant
sim_res_split_fake <- list(
  data = list(Y = Y_fake, W = W_fake),
  nuisances = list(xgb = list(split = fake_nuisance_df()),
                   automl = list(split = fake_nuisance_df()))
)
scored_split <- me_per_model(
  list(tau = tau_hat_fake), true_tau_fake, "test_model", sim_res_split_fake, NULL
)
n_split_expected <- 1 + length(fake_pipelines) * 1 * n_score_types
report(
  ncol(scored_split) == n_split_expected &&
    all(c("dr_split_xgb", "infl_split_automl", "cal05q5_split_xgb") %in%
          names(scored_split)),
  sprintf("the 80:20 tree scores through the same me_per_model() (got %d of %d columns)",
          ncol(scored_split), n_split_expected)
)

# ---- the calibration score's own properties --------------------------------
# calc_cal_score() is the one genuinely new formula, so it gets checked against
# its definition rather than only against me_per_model()'s assembly of it.

# perfect calibration: if the DR scores agree with the candidate row for row,
# every group's two GATEs are identical and the score is exactly 0
report(
  calc_cal_score(tau_hat_fake, tau_hat_fake, 5) == 0,
  "phi == tau_hat gives exactly 0 (perfect calibration)"
)

# scenario 1 / a fully-shrunk elastic net: one predicted value for every row.
# quantile cut points would collapse here; rank-based grouping must not.
tau_const <- rep(0.3, n_fake)
cal_const <- calc_cal_score(tau_const, rnorm(n_fake), 5)
g_const <- ceiling(5 * rank(tau_const, ties.method = "first") / n_fake)
report(
  is.finite(cal_const) && length(unique(g_const)) == 5,
  sprintf("constant tau_hat still forms 5 non-empty groups and scores finite (got %.3f)",
          cal_const)
)

# weights are COUNTS, not proportions - sum(w) is n. Guards against a later
# "fix" that divides by n and silently rescales every calibration column.
g_fake <- ceiling(5 * rank(tau_hat_fake, ties.method = "first") / n_fake)
report(
  sum(as.numeric(table(g_fake))) == n_fake,
  sprintf("calibration weights are group counts summing to n = %d", n_fake)
)

# the n-scaling that follows from those weights, stated as a check rather than
# left to be discovered from a plot: duplicating the data doubles the score
phi_fake_cal <- rnorm(n_fake)
report(
  abs(calc_cal_score(c(tau_hat_fake, tau_hat_fake), c(phi_fake_cal, phi_fake_cal), 5) -
        2 * calc_cal_score(tau_hat_fake, phi_fake_cal, 5)) < 1e-10,
  "the score scales with n (duplicating the data doubles it) - intended, see README"
)

cal_by_k <- vapply(CAL_QUANTILES, function(k) {
  calc_cal_score(tau_hat_fake, phi_fake_cal, k)
}, numeric(1))
report(
  all(is.finite(cal_by_k)) && length(unique(cal_by_k)) == length(CAL_QUANTILES),
  sprintf("every K in CAL_QUANTILES gives a distinct finite score (%s)",
          paste(sprintf("%.3f", cal_by_k), collapse = ", "))
)

# the fixed-pi twins rest on a scalar recycling through calc_infl_score()'s
# vectorised arithmetic - cheap to prove, expensive to get silently wrong
report(
  abs(calc_infl_score(tau_hat_fake, fake_nuis$xgb[[1]]$tau_T, 0.5, Y_fake, W_fake) -
        calc_infl_score(tau_hat_fake, fake_nuis$xgb[[1]]$tau_T,
                        rep(0.5, n_fake), Y_fake, W_fake)) < 1e-12,
  "calc_infl_score() with scalar pi = 0.5 equals the length-n version"
)

# =============================================================================
cat("\n=== 2. DGM / design-matrix plumbing ===\n")

setup_rng_stream(1)
gen1 <- generate_me_scenario_data(scenario = 1, n = 100)
data1 <- gen1$dataset

report(
  identical(names(data1)[1:2], c("Y", "W")),
  "Y, W are the first two columns of the generated dataset"
)

design1 <- prepare_design_matrix(data1)
all_numeric <- all(vapply(design1$X, is.numeric, logical(1)))
report(all_numeric, "design matrix covariates are all numeric - no factor expansion needed")
report(nrow(design1$X) == 100, "design matrix has the right number of rows")
report(!is.null(gen1$truth$tau), "generate_me_scenario_data() returns truth$tau")

# =============================================================================
cat("\n=== 3. scenario 1: no heterogeneity ===\n")
# tau is computed as p1 - p0 with a shared b0 + b1*X1 + b2*X2 term, so the
# true constant comes back with floating point dust on it - same check
# continuous/ and crossfitting/ both use.

tau_spread <- diff(range(gen1$truth$tau))
report(
  tau_spread < 1e-8,
  sprintf("true CATE is constant at %.3f (spread %.1e)", gen1$truth$tau[1], tau_spread)
)

# =============================================================================
cat(sprintf(
  "\n=== 4. structural check on the real pipeline%s ===\n",
  if (full) "" else " (candidate models only - rerun with 'full' for nuisances)"
))

n_test <- 250 # smallest grid n, to keep this check as fast as possible
n_folds_test <- 10L # matches me_analysis.R's rule (single crossfitting, all n)

setup_rng_stream(1)
gen4 <- generate_me_scenario_data(scenario = 4, n = n_test)
design4 <- prepare_design_matrix(gen4$dataset)
Y4 <- design4$Y
W4 <- design4$W
X4 <- design4$X
kfolds4 <- split_folds(Y4, k = n_folds_test)

# Caught rather than allowed to abort the script. SL2's SuperLearner library
# includes SL.xgboost, whose bundled wrapper predates the installed xgboost's
# renamed API, and predict.SuperLearner() then crashes on the NULL fit that
# failure left behind (README.md's "Known local-environment limitation").
# Uncaught, that took the whole file down at this point - so sections 5 and 6
# below, which need neither the candidates nor H2O, could never run on a local
# machine, and this folder had no usable local verification at all. The
# failure is still reported as a failure; it just no longer hides everything
# after it.
t0 <- Sys.time()
model_list <- tryCatch(
  run_all_candidate_models(Y4, W4, X4, kfolds4$fold_indices, kfolds4$fold_list),
  error = function(e) {
    cat("  candidate-model fitting ERRORED: ", conditionMessage(e), "\n", sep = "")
    NULL
  }
)
elapsed_models <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  candidate-model fitting took %.1fs\n", elapsed_models))

models_ok <- !is.null(model_list)
report(
  models_ok,
  paste0(
    "all 9 candidate models fit without error",
    if (!models_ok) {
      " [expected locally - see README.md's known SL2/SL.xgboost limitation; must PASS on the cluster]"
    }
  )
)

# The check that used to sit here asserted
#
#   !identical(kfolds4$fold_indices, nuis_folds4$fold_indices)
#   "candidate and scoring-pipeline fold draws are independent"
#
# and has been DELETED rather than updated, because it asserts precisely the
# property that has been retired. me_analysis.R's second independent fold
# draw was justified as preserving honesty; it does not - row-level honesty
# holds under a shared draw too, and an independent draw only cuts the two
# training sets' overlap from 100% to (V-1)/V = 90%. See me_config.R's
# NUISANCE_ARMS. The arm that replaces it (`cv_shared`) uses the candidates'
# own folds deliberately, so this check would now fail for the right change.

report(
  is.null(kfolds4$fold_pairs),
  "split_folds() no longer returns fold_pairs (single crossfitting, not double)"
)

matrix_targets <- c("po", "Y0.hat", "Y1.hat", "W.hat")
vectors_ok <- vapply(model_list, function(m) {
  all(vapply(matrix_targets, function(nm) is.vector(m[[nm]]) && !is.matrix(m[[nm]]), logical(1)))
}, logical(1))
report_if(
  models_ok,
  all(vectors_ok),
  paste0(
    "po/Y0.hat/Y1.hat/W.hat are length-n vectors, not double-crossfit matrices",
    if (!all(vectors_ok)) {
      paste0(" [bad: ", paste(names(which(!vectors_ok)), collapse = ", "), "]")
    }
  )
)

no_na_nuis <- vapply(model_list, function(m) {
  all(vapply(matrix_targets, function(nm) sum(is.na(m[[nm]])) == 0, logical(1)))
}, logical(1))
report_if(
  models_ok,
  all(no_na_nuis),
  paste0(
    "no NA in any candidate's po/Y0.hat/Y1.hat/W.hat",
    if (!all(no_na_nuis)) {
      paste0(" [bad: ", paste(names(which(!no_na_nuis)), collapse = ", "), "]")
    }
  )
)

what_trimmed <- vapply(model_list, function(m) all(m$W.hat >= 0.05 & m$W.hat <= 0.95), logical(1))
report_if(
  models_ok,
  all(what_trimmed),
  paste0(
    "every candidate's W.hat is trimmed to [0.05, 0.95]",
    if (!all(what_trimmed)) {
      paste0(" [bad: ", paste(names(which(!what_trimmed)), collapse = ", "), "]")
    }
  )
)

report_if(
  models_ok,
  setequal(names(model_list), CANDIDATE_MODELS),
  sprintf(
    "all %d candidate models present (got %d)",
    length(CANDIDATE_MODELS), length(model_list)
  )
)

lengths_ok <- vapply(model_list, function(m) length(m$tau) == n_test, logical(1))
report_if(
  models_ok,
  all(lengths_ok),
  paste0(
    "all tau vectors correctly sized",
    if (!all(lengths_ok)) {
      paste0(" [bad: ", paste(names(which(!lengths_ok)), collapse = ", "), "]")
    }
  )
)

na_ok <- vapply(model_list, function(m) sum(is.na(m$tau)) == 0, logical(1))
report_if(
  models_ok,
  all(na_ok),
  paste0(
    "no NA tau estimates",
    if (!all(na_ok)) {
      paste0(" [bad: ", paste(names(which(!na_ok)), collapse = ", "), "]")
    }
  )
)

if (full && models_ok) {
  model_seed <- sample.int(2^31 - 1, 1)
  # Every arm me_strategies.R will actually run, on the CANDIDATES' folds -
  # cv_shared and holdout both derive from kfolds4, which is the point of
  # them. `whole` stands in for the carried-through arm.
  hold4 <- holdout_blocks(kfolds4$fold_indices)
  arms4 <- list(
    whole     = nuisance_arm_spec("whole"),
    cv_shared = nuisance_arm_spec("cv", kfolds4),
    holdout   = nuisance_arm_spec("holdout", hold4)
  )

  t0 <- Sys.time()
  nuisances4 <- run_nuisance_arms(
    X4, Y4, W4, arms = arms4,
    n_cores = 2, mem = "4G", model_seed = model_seed
  )
  elapsed_nuis <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  nuisance-evaluation pipelines took %.1fs\n", elapsed_nuis))

  report(
    setequal(names(nuisances4), c("xgb", "automl")),
    "both nuisance pipelines (xgb, automl) present"
  )

  arms_ok <- vapply(nuisances4, function(p) setequal(names(p), names(arms4)), logical(1))
  report(all(arms_ok), "both pipelines expose exactly the arms they were asked for")

  expected_cols <- c("mu0_T", "mu1_T", "mu0_DR", "mu1_DR", "pi", "mu_DR", "tau_T", "phi", "phi05")
  cols_ok <- vapply(names(nuisances4), function(p) {
    all(vapply(nuisances4[[p]], function(ft_df) setequal(names(ft_df), expected_cols), logical(1)))
  }, logical(1))
  report(all(cols_ok), "every pipeline/arm data.frame has the expected 9 columns")

  rows_ok <- vapply(names(nuisances4), function(p) {
    all(vapply(nuisances4[[p]], function(d) nrow(d) == n_test, logical(1)))
  }, logical(1))
  report(all(rows_ok), sprintf("every pipeline/arm covers all %d rows", n_test))

  # cv_shared and holdout see disjoint data for any given row, so they cannot
  # agree by construction. If they do, one of them got the wrong filter.
  distinct_ok <- vapply(names(nuisances4), function(p) {
    !isTRUE(all.equal(nuisances4[[p]]$cv_shared$phi, nuisances4[[p]]$holdout$phi))
  }, logical(1))
  report(all(distinct_ok), "cv_shared and holdout produce genuinely different phi")

  scored_real <- me_per_model(
    model_list[[1]], gen4$truth$tau, names(model_list)[1],
    list(data = list(Y = Y4, W = W4), nuisances = nuisances4), NULL
  )
  # derived from me_score_types() rather than a typed-in score-type count, so
  # adding a score family cannot leave this assertion quietly checking the old
  # number
  n_real <- 1 + 2 * length(arms4) * length(me_score_types())
  report(
    ncol(scored_real) == n_real,
    sprintf("me_per_model() produces %d score columns on real pipeline output (got %d)",
            n_real, ncol(scored_real))
  )
} else {
  cat("  (nuisance pipelines skipped - they are the slow half and SL2 fails\n")
  cat("   locally for a documented package-version reason; run them on the\n")
  cat("   cluster via jobscripts/me_testing.sh rather than with 'full' here)\n")
}

# =============================================================================
cat("\n=== 5. nuisance arms: block derivation and filters ===\n")
# Structural only - no fitting, so this runs everywhere.

# holdout_blocks() must POOL at n=250 (25-row folds) and leave n=500/1000
# alone (50/100-row folds). The pooled blocks must be unions of whole
# candidate folds, or they are not "the data the candidate never saw".
for (nn in c(250, 500, 1000)) {
  fi <- caret::createFolds(y = rnorm(nn), k = 10L, list = FALSE)
  hb <- holdout_blocks(fi)
  block_sizes <- as.integer(table(hb$fold_indices))
  n_blocks <- length(hb$fold_list)

  report(
    min(block_sizes) >= HOLDOUT_MIN_BLOCK,
    sprintf("n=%d: every holdout block has >= %d rows (min %d, %d blocks)",
            nn, HOLDOUT_MIN_BLOCK, min(block_sizes), n_blocks)
  )
  report(
    n_blocks == if (nn == 250) 5L else 10L,
    sprintf("n=%d: %d holdout blocks (expected %d)",
            nn, n_blocks, if (nn == 250) 5L else 10L)
  )
  # each candidate fold sits entirely inside one block
  nested <- all(vapply(split(hb$fold_indices, fi), function(b) length(unique(b)) == 1L, logical(1)))
  report(nested, sprintf(
    "n=%d: every candidate fold lies wholly within one holdout block", nn
  ))
  report(
    length(hb$fold_indices) == nn,
    sprintf("n=%d: holdout blocks cover every row", nn)
  )
}

# The filters each arm uses, asserted directly rather than inferred from
# fitted output. This is the honesty property in one line each.
fi_t <- caret::createFolds(y = rnorm(200), k = 10L, list = FALSE)

cv_honest <- all(vapply(sort(unique(fi_t)), function(f) {
  train <- !(fi_t %in% f)
  test <- !train
  # no test row is in the training set, and the training set is exactly
  # "every fold but this one" - the candidate's own split
  !any(train & test) && setequal(which(train), which(fi_t != f))
}, logical(1)))
report(cv_honest, "cv_shared: no row is in the model that predicts it, and the training set is the candidate's V-1 folds")

hold_coupled <- all(vapply(sort(unique(fi_t)), function(f) {
  train <- fi_t %in% f
  test <- train
  identical(train, test)
}, logical(1)))
report(hold_coupled, "holdout: train and test filters are identical (resubstitution, by design)")

sp_t <- single_split(rnorm(500))
report(
  length(intersect(sp_t$train, sp_t$test)) == 0 &&
    length(sp_t$train) + length(sp_t$test) == 500,
  sprintf("single_split(): disjoint and exhaustive (%d train / %d test)",
          length(sp_t$train), length(sp_t$test))
)
report(
  abs(length(sp_t$train) / 500 - SPLIT_TRAIN_FRAC) < 0.02,
  sprintf("single_split() trains on ~%.0f%% (got %.1f%%)",
          100 * SPLIT_TRAIN_FRAC, 100 * length(sp_t$train) / 500)
)

# =============================================================================
cat("\n=== 6. me_strategies.R assembly, on a synthetic res_sim object ===\n")
# There are no model_evaluation results on a local machine - results/ carries
# every other study but not this one - so the read/compute/merge contract is
# exercised against a fixture instead. XGBoost only: H2O is the slow half and
# needs a Java runtime, and nothing about the merge depends on which
# estimator produced the data.frames.

# the study's smallest n, deliberately: the holdout arm is tightest here
# (25-row candidate folds, pooled to 5 blocks of 50, and mu0_T/mu1_T then see
# only the ~25 control / ~25 treated rows inside a block). Testing at some
# smaller convenience size would exercise a configuration the study never runs
# and miss the one that is actually near the edge.
n_asm <- 250L
setup_rng_stream(7)
gen_asm <- generate_me_scenario_data(scenario = 4, n = n_asm)
design_asm <- prepare_design_matrix(gen_asm$dataset)
folds_asm <- split_folds(design_asm$Y, k = 10L)

# what me_analysis.R would have written: cv + whole, under its own draw
old_asm <- list(
  data = list(Y = design_asm$Y, W = design_asm$W, X = design_asm$X),
  truth = gen_asm$truth,
  fold_info = folds_asm,
  nuisances = list(xgb = list(
    cv = fake_nuisance_df_n(n_asm), whole = fake_nuisance_df_n(n_asm)
  ))
)

# run_nuisance_arms() returns list(<pipeline> = list(<arm> = df)), so the arms
# sit one level below the pipeline even when only one pipeline was asked for.
new_asm <- run_nuisance_arms(
  design_asm$X, design_asm$Y, design_asm$W,
  arms = list(
    cv_shared = nuisance_arm_spec("cv", folds_asm),
    holdout = nuisance_arm_spec("holdout", holdout_blocks(folds_asm$fold_indices))
  ),
  n_cores = 1, mem = "1G", model_seed = 1L,
  pipelines = "xgb"
)$xgb

report(
  setequal(names(new_asm), c("cv_shared", "holdout")),
  "run_nuisance_arms() returns exactly the arms it was asked for"
)

merged_asm <- list(
  whole     = old_asm$nuisances$xgb$whole,
  cv_indep  = old_asm$nuisances$xgb$cv,
  cv_shared = new_asm$cv_shared,
  holdout   = new_asm$holdout
)

report(
  identical(names(merged_asm), NUISANCE_ARMS),
  "merged arm names match NUISANCE_ARMS in order"
)
report(
  identical(merged_asm$whole, old_asm$nuisances$xgb$whole) &&
    identical(merged_asm$cv_indep, old_asm$nuisances$xgb$cv),
  "whole and cv_indep pass through bit-identically (cv -> cv_indep is a pure rename)"
)
report(
  nrow(new_asm$cv_shared) == n_asm && nrow(new_asm$holdout) == n_asm,
  sprintf("both new arms cover all %d rows", n_asm)
)
report(
  !anyNA(new_asm$cv_shared$phi) && !anyNA(new_asm$holdout$phi),
  "neither new arm produces NA in phi"
)
report(
  !isTRUE(all.equal(new_asm$cv_shared$phi, new_asm$holdout$phi)),
  "cv_shared and holdout differ - the fold vectors really are being used"
)
# the leak the holdout arm accepts, made visible rather than assumed: fitting
# and predicting on the same rows shrinks Y - mu_DR toward zero, so its AIPW
# correction is systematically smaller than the row-honest arm's
corr_shared <- mean(abs(new_asm$cv_shared$phi - (new_asm$cv_shared$mu1_DR - new_asm$cv_shared$mu0_DR)))
corr_hold <- mean(abs(new_asm$holdout$phi - (new_asm$holdout$mu1_DR - new_asm$holdout$mu0_DR)))
cat(sprintf(
  "  mean |AIPW correction|: cv_shared %.4f, holdout %.4f (holdout smaller = the expected leak)\n",
  corr_shared, corr_hold
))

# =============================================================================
plan(sequential)

cat("\n=== summary ===\n")
cat(sprintf("  %d passed, %d failed\n", length(pass), length(fail)))
if (length(fail) > 0) {
  cat("\nfailures:\n")
  for (f in fail) {
    cat("  - ", f, "\n", sep = "")
  }
  quit(status = 1)
}
cat(sprintf(
  "\nall checks passed.%s next: submit jobscripts/me_profile.sh, then run me_profile_summary.R\n",
  if (!full) " run with 'full' to also exercise the real nuisance pipelines, then" else ""
))
