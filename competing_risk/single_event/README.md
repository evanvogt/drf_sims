# Single-event survival

CATE estimation on the restricted mean survival time with one event and no
competing risk. The event is `competing_risk/`'s event 1 with event 2 removed,
so this study's estimand is that study's net RMST1 (`tau_RMST1_cs`), and the
two can be read side by side. Full design in `ADEMP.md`.

| | |
|---|---|
| scenarios | 1 null, 2 constant (hazard scale), 3 heterogeneous |
| n | 500 |
| censoring | TRUE (uniform + administrative at 180), FALSE (administrative only) |
| runs | 500 |
| array | **3,000 jobs** (`se_1.sh`) |
| horizon | 28 |
| results | `../results/single_event/scenario_<k>/500/censor_<TRUE\|FALSE>/` |

Results sit in `results/single_event/`, not under `results/competing_risk/`,
because the parent README's archive step tars and removes that whole tree.

On the cluster, from `competing_risk/single_event/jobscripts/` (run
`Rscript make_log_dirs.R` from the repo root first, for `logs_1/` and
`logs_rerun/`):

```bash
qsub se_1.sh                # 1-3000
Rscript ../se_check.R       # writes failed_ids.txt, points se_rerun.sh at it
qsub se_rerun.sh            # only if the check found failures
qsub se_collect.sh
qsub se_metrics.sh
```

Local smoke test, from `competing_risk/single_event/`:
`Rscript se_analysis.R 3` (scenario 3, censoring on, run 1), then
`Rscript se_analysis.R 6` (the same run, censoring off). Each takes about 1.5
minutes locally with the default 2 workers. Don't run `se_check.R` locally: it
would rewrite `failed_ids.txt` and `se_rerun.sh` for every run you haven't
produced.

`Rscript se_dgm_check.R` (about 30 s) checks that the truth matches the
closed form and the parent's `tau_RMST1_cs`, and that the pairing below holds.
It also prints the event mix and population truths quoted in `ADEMP.md`.

## Arms

| learner | random forest on pseudo-values | SuperLearner on pseudo-values | RSF on (Y, D) |
|---|---|---|---|
| DR-learner | `pseudo_dr_whole_oob` | `sl_dr_whole` | `rsf_dr_oob`, `rsf_dr_scf` |
| T-learner | `pseudo_t_whole_oob` | `sl_t_whole` | `rsf_t_oob`, `rsf_t_scf` |

Plus `csf` (causal survival forest) and `pseudo_cf_whole_oob` (causal forest on
pseudo-values).

- **The T-learners cost nothing extra.** `t_from_dr()` reads μ̂₁ − μ̂₀ off the
  DR-learner's own per-arm outcome models, so each T-learner differs from its
  DR-learner in the second stage alone.
- **Reused, not copied.** Everything except the RSF outcome model,
  `pseudo_rmst()` and `t_from_dr()` comes unchanged from
  `competing_risk/surv_models.R`, so a change there changes this study too.
  `R/regression_check.R` has a `single_event` entry for that.
- **The RSF outcome model is a survival-family forest.** The parent's
  `rsf_fit`/`rsf_rmtl` would take their single-cause warning branch on every
  fit, so `se_models.R` has its own (`rsf_se_fit`, `rsf_rmst`), on the same
  pattern.

## Things to expect when reading results

- **Without censoring the pseudo-values are exactly min(Y, 28).** The
  pseudo-value arms are then regressions on the observed truncated time, which
  the smoke test confirms to 1e-12.
- **Censoring is light.** About 10% of each arm is censored before 28, so the
  censored and uncensored cells may not differ much.
- **Scenario 2's correlation and C-statistic mean little.** The truth's SD
  there is 0.10, all from X1 and X2, so read bias and RMSE.
- **Scenario 1 is a true null**, so `cate_metrics()`'s hard-coded scenario-1
  convention (correlation 0, C 0.5) is correct here. `tau_sd` measures the
  spurious heterogeneity each arm reports.
- **"SuperLearner failed for W.hat. Using mean(W)" warnings are expected.** In
  an RCT the propensity library often gets zero weight everywhere. The
  failsafe is the parent's, and mean(W) is the right answer anyway.

## Pairing with `competing_risk/`

A run here draws W, X and U exactly as `competing_risk/`'s ρ = 0 run with the
same run index does. Only the parent's `cause` draw is missing. So:

- the two studies' datasets share treatment and covariates unit by unit;
- each unit's event time here is at least its all-cause time there.

So per-run differences between the studies (this study's `csf` against the
parent's `csf_cs` Event 1, for example) are paired. C is not shared: the
parent draws `cause` before C.

Since 2026-10-06 the parent's primary analysis is ρ = 0.5 and ρ = 0 is its
sensitivity analysis, so this study pairs with the parent's sensitivity runs.
That was a deliberate choice, to avoid rerunning this study at ρ = 0.5.
Against the parent's primary (ρ = 0.5) results a comparison is side by side,
not paired.

## Files

| file | |
|---|---|
| `se_config.R` | the grid |
| `se_dgm.R` | DGM and closed-form truth |
| `se_models.R` | `all_cate_se_models()`, plus the single-event pieces |
| `se_analysis.R` | one array index |
| `se_check.R`, `se_collect.R`, `se_metrics.R` | the pipeline |
| `se_dgm_check.R` | DGM assertions and the ADEMP numbers |
| `se_results.qmd` | the results report (render once the results are back) |
