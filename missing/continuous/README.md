# Missing covariates — continuous outcome

The main missing-data design. See `missing/README.md` for the mechanisms,
handling methods and the shared bug fixes (B, C, D — all of which originated
here).

| | |
|---|---|
| array | **9,900 jobs**, split across `cts_miss_1.sh` (1–5000) and `cts_miss_2.sh` (5001–9900) |
| results | `../results/missing/continuous/scenario_<k>/<n>/<type>/<prop>/<mechanism>/<method>/` |
| metrics | `cts_miss_metrics.RDS`, including `rel_efficiency` against the `complete_data` arm; plus `cts_miss_true_cate_tests.RDS`, the true-CATE HTE test evaluation |

## Files

| file | role |
|---|---|
| `cts_miss_config.R` | the grid — **the** definition, never filtered |
| `cts_miss_dgms.R` | names the `continuous_missing` scenario set + the missingness machinery |
| `cts_miss_models.R` | `family = gaussian()`, `profile = "missing"`, `ipw` threaded through |
| `cts_miss_analysis.R` | array entry point |
| `cts_miss_patch.R` | one-off: back-fills the `dr_random_forest` HTE tests into finished results (`R/patch_hte_tests.R`). Array entry point for `jobscripts/cts_miss_patch.sh` |
| `cts_miss_patch_check.R` | audits that back-fill — which combinations it failed on and why. See `missing/binary/README.md`, where the failure it was written for happened; this study's patch is complete, so a run here should report nothing to re-run |
| `cts_miss_results.R` / `.qmd` | every metric, to `results/all_figures/` — the diagnostic counterpart to the chapter script `results_processing/thesis_figures/miss_cts.R` |
| `mi_scratch.R` | exploratory, unmaintained |
| `results_cts_miss.{R,Rmd,html}` | **superseded** by `cts_miss_results.R`/`.qmd` — kept for reference, renamed to the repo's legacy convention (cf. `continuous/results_cts.R`). Points at the old HPC path `/rds/general/...` and at `results/new_format/metrics_cts_miss_df.RDS`, neither of which exists any more; scenarios are `"scenario_5"` strings, and only bias and MSE are plotted. It was called `cts_miss_results.*`, so its `.html` was being clobbered by the new `.qmd`'s render output. |

The `ipw` argument is the only thing separating this study's estimators from
`continuous/`: it becomes grf `sample.weights` and SuperLearner `obsWeights`.
`ipw = NULL` reproduces the unweighted path exactly.

## Gotchas

**Scenario numbers are not the main study's.** Scenario `k` here maps to
main-study scenario 1, 2, 4, 8, 9.

**`multiple_imputation` returns a list of 50 datasets, not one.** The analysis
script fits each and Rubin-combines with `combine_mi()`; only
`causal_forest`, `dr_random_forest` and `dr_semi_oracle` are combined.

**The row-dropping methods drop rows from the truth too.** `complete_cases` and
`IPW` return `retained_indices`, and `generate_and_process_data()` subsets
`truth` to match — otherwise estimates and truth would be misaligned.

**`cts_miss_true_cate_tests.RDS` runs the BLP and independence tests on the
true CATE and true nuisances instead of an estimator's** (`truth$tau`,
`truth$p0`, `W.hat = 0.5` — see `continuous/README.md`'s "True-CATE HTE test
evaluation"), one row per (scenario, n, type, prop, mechanism, method, run),
no per-model dimension. `method == "multiple_imputation"` rows are `NA`/`NA`
there, same reason as the estimated-CATE gap below: `data` is a list of 50
imputed data.frames, with no single covariate matrix to test against.

**`dr_random_forest` used to carry no BLP or independence test in this study**,
unlike `continuous/`, so `BLP_p` was `NA` for that one model. Copy-paste drift
rather than a decision, and the decision has been taken: every model carries the
tests where possible. `PROFILES$missing` (`R/cate_models.R`) now sets
`dr_rf_tests = TRUE`, and the finished results were back-filled in place by
`R/patch_hte_tests.R` rather than re-run — the tests are deterministic and every
input survives in the saved files. Written up once, in
`missing/binary/README.md`, since `profile = "missing"` covers both studies;
that write-up also covers the one field the patch cannot recover
(`dr_random_forest$variance`, which nothing reads).

Patch this study only after `cts_miss_rerun.sh` has cleared the outstanding runs
and `cts_miss_check.R` reports 9,900/9,900, so it is one clean pass:

```bash
Rscript cts_miss_patch.R dry                       # report only, writes nothing
qsub    jobscripts/cts_miss_patch.sh               # 1-99, one combination each
Rscript cts_miss_patch_check.R                     # did every element land?
```

That last step is not ceremony: the binary study's patch job lost ten of its 99
array elements silently, because a manifest is only written once a combination
finishes and `check_all.R` counts manifest rows. See "Checking the back-fill
landed" in `missing/binary/README.md`.

## Running it

```bash
qsub missing/continuous/jobscripts/cts_miss_1.sh   # 1-5000
qsub missing/continuous/jobscripts/cts_miss_2.sh   # 5001-9900
Rscript missing/continuous/cts_miss_check.R
qsub missing/continuous/jobscripts/cts_miss_patch.sh    # 1-99, the HTE back-fill
qsub missing/continuous/jobscripts/cts_miss_collect.sh
qsub missing/continuous/jobscripts/cts_miss_metrics.sh
```

The patch step goes **before** collect: collect reads the per-run files into
`cts_miss_all.RDS`, so patching afterwards would leave the collected copy
carrying the old, testless `dr_random_forest`. It is a one-off — once these
results are patched and `PROFILES$missing` is set, future runs need only the
usual steps.

Then, for the figures:

```bash
Rscript missing/continuous/cts_miss_results.R           # every metric, to results/all_figures/
quarto render missing/continuous/cts_miss_results.qmd   # the same, as a browsable report
```

To run only the `complete_data` reference arm, take
`grid_indices(study, method = "complete_data")` (1,100 indices) and submit those
— do **not** filter the grid.

## Status

**Re-runs required** — for the crossfitting strategy change to
`R/cate_models.R` (see root README Methods/Status), which moves all five
estimator arms, and separately for bug F (`dr_superlearner` only). The DGM is
unaffected by the *bug ledger*. Bugs B and C were collection/metrics
problems, so re-running `cts_miss_collect.R` and `cts_miss_metrics.R` over
the existing per-run files recovers the mechanisms that were being missed
and populates `rel_efficiency` — no cluster time needed for those two.
