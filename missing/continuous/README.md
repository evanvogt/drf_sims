# Missing covariates — continuous outcome

The main missing-data design. See `missing/README.md` for the mechanisms,
handling methods and the shared bug fixes (B, C, D — all of which originated
here).

| | |
|---|---|
| array | **9,900 jobs**, split across `cts_miss_1.sh` (1–5000) and `cts_miss_2.sh` (5001–9900) |
| results | `../results/missing/continuous/scenario_<k>/<n>/<type>/<prop>/<mechanism>/<method>/` |
| metrics | `cts_miss_metrics.RDS`, including `rel_efficiency` against the `complete_data` arm |

## Files

| file | role |
|---|---|
| `cts_miss_config.R` | the grid — **the** definition, never filtered |
| `cts_miss_dgms.R` | names the `continuous_missing` scenario set + the missingness machinery |
| `cts_miss_models.R` | `family = gaussian()`, `profile = "missing"`, `ipw` threaded through |
| `cts_miss_analysis.R` | array entry point |
| `mi_scratch.R` | exploratory, unmaintained |

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

**`dr_random_forest` carries no BLP or independence test in this study**, unlike
`continuous/`, so `BLP_p` is `NA` for that one model. This looks like copy-paste
drift rather than a decision; it is preserved (see the profile table in
`R/README.md`) but is worth revisiting.

## Running it

```bash
qsub missing/continuous/jobscripts/cts_miss_1.sh   # 1-5000
qsub missing/continuous/jobscripts/cts_miss_2.sh   # 5001-9900
Rscript missing/continuous/cts_miss_check.R
qsub missing/continuous/jobscripts/cts_miss_collect.sh
qsub missing/continuous/jobscripts/cts_miss_metrics.sh
```

To run only the `complete_data` reference arm, take
`grid_indices(study, method = "complete_data")` (1,100 indices) and submit those
— do **not** filter the grid.

## Status

**Re-runs for bug F only** (`dr_superlearner`). The DGM is unaffected by the bug
ledger. Bugs B and C were collection/metrics problems, so re-running
`cts_miss_collect.R` and `cts_miss_metrics.R` over the existing per-run files
recovers the mechanisms that were being missed and populates `rel_efficiency` —
no cluster time needed for those two.
