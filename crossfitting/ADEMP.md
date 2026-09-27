# ADEMP — crossfitting comparison

## Aims

Compare double crossfitting with cheaper sample-splitting schemes (single
crossfit, out-of-bag) for the nuisance and second-stage models of DR-learner
and causal forest CATE estimators: accuracy, generalisation and run time.

A pilot (`confidence_intervals/`) also compares interval estimates across
the RF and causal forest arms.

## Data-generating mechanisms

The `continuous/` DGM (see `continuous/ADEMP.md`), at:

| | |
|---|---|
| scenarios | 1 (no HTE), 4 (cosine), 6 (two variables), 8 (single effects + interaction) |
| n | 500 (training) + independent test sample of 2000 from the same DGM |
| repetitions | 100 per scenario |

All arms in a replicate share one fold assignment (V = 10).

CI pilot: scenarios 1, 4, 8; n = 500; 50 repetitions.

## Estimands

- Unit-level CATE on the training sample.
- CATE at the covariates of the independent test sample.

## Methods

DR-learner, random forest (7 arms):

| arm | nuisances | second stage |
|---|---|---|
| `dcf` | double crossfit over fold pairs | crossfit, same folds |
| `scf_scf` | single crossfit | crossfit, same folds |
| `scf_oob` | single crossfit | whole sample, OOB |
| `scf_oob_t` | single crossfit, T-learner | whole sample, OOB |
| `oob_oob` | whole sample, OOB, T-learner | whole sample, OOB |
| `oob_oob_s` | whole sample, OOB, S-learner (`X.orig` workaround) | whole sample, OOB |
| `oob_oob_manual` | whole sample, OOB, S-learner (manual tree loop) | whole sample, OOB |

DR-learner, SuperLearner: `dcf`, `scf_scf`.

Causal forest: `cf_dcf` (double-crossfit `Y.hat`/`W.hat`, fold-wise forest),
`cf_scf` (single-crossfit nuisances, fold-wise), `cf_full_oob` (single-crossfit
nuisances, whole-sample OOB), `cf_default` (plain grf `causal_forest`).

Propensities trimmed to [0.05, 0.95] in every arm.

CI pilot: the 11 RF / causal forest arms, each with a half-sample bootstrap
simultaneous band (`CI_boot = 200`, `CI_sf = 0.5`, α = 0.05). OOB arms carry
two bootstrap variants (`half_boot`, `half_boot_out`) and grf's pointwise
normal interval (`grf_normal`).

## Performance measures

Per arm, on both the training and test sample (`R/metrics.R::cate_metrics`):

- bias, ATE bias, relative biases
- MSE, RMSE, MAE; `mse_test_single` (each fold model scored separately on the test set, then averaged)
- Pearson and Spearman correlation, sign accuracy
- run time: nuisance, second stage, total

CI pilot: marginal coverage, simultaneous coverage (nominal 0.95), mean interval length.
