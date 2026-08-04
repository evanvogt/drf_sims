# Missing covariates

What happens to CATE estimation when covariates are missing, and which handling
method costs least?

Three studies:

| folder | question |
|---|---|
| `continuous/` | the main design — continuous outcome, 8 handling methods × 3 mechanisms |
| `binary/` | the same design on a logit scale |
| `ci_example/` | how to build a *confidence interval* after multiple imputation |

## The design (continuous and binary)

| | |
|---|---|
| scenarios | 1, 2, 4, 5 — a **renumbered** five-scenario set, not a subset of the main ten |
| n | 500 |
| type | `both` (prognostic and predictive covariates amputated) |
| prop | 0.3 |
| mechanism | MAR, MNAR, MNAR-Y |
| method | 8 handling methods + `complete_data` reference |
| runs | 100 |
| array | **9,900 jobs** (scenario 1 has no MNAR-Y, so 100 rows are dropped) |

Scenario numbering here is **not** the main study's. Scenario `k` corresponds to
main-study scenario 1, 2, 4, 8, 9 respectively — so "scenario 4" means something
different in this folder than in `continuous/`.

### Mechanisms

- **MAR** — missingness depends on the observed covariates
- **MNAR** — driven entirely by an unobserved `U`
- **MNAR-Y** — `U` also enters the treatment effect, so missingness is related to
  the outcome. Not defined for scenario 1, which has no treatment effect
  heterogeneity to relate to.

`missing/ci_example/` still calls these `AUX` / `AUX-Y`. Both spellings are
accepted and normalised in `R/missingness.R`.

### Handling methods

`complete_cases`, `mean_imputation`, `missforest`, `regression`,
`missing_indicator`, `IPW`, `multiple_imputation`, `none` (let the estimator
handle it), plus `complete_data` — a reference arm with **no missingness
introduced at all**, so the others can be scored against complete-data
performance via `rel_efficiency`.

## Bugs fixed here

**Bug B — `cts_miss_collect.R` looked for directories that could not exist.** It
used `mechanism = c("MAR", "AUX", "AUX-Y")` while the DGM, analysis and check
scripts all used `MNAR`/`MNAR-Y`. Two of the three mechanisms were silently
collected as empty. `missing/binary/` was already consistent.

**Bug C — `rel_efficiency` was `NA` everywhere.** The metrics scripts built the
reference from `method == "complete_data"`, which the collect grid did not
contain, so the join found nothing.

**Bug D — the array index meant two different things.** Both analysis scripts
appended `filter(method == "complete_data")` *after* `expand.grid`, renumbering
every row. A `failed_ids.txt` written by the check script would have resubmitted
the wrong parameters. The grid now lives in `<prefix>_config.R` and is never
filtered; to run one arm, select indices instead:

```r
idx <- grid_indices(study, method = "complete_data")
```

This is also why the jobscript ranges changed: `1-1100` was the *post-filter*
size, and `cts_miss_2.sh`'s `10001-12200` referred to indices that never existed.

## Status

`missing/continuous` — re-runs for bug F (`dr_superlearner` only).

`missing/binary` — re-runs entirely; see its own README, it had three further
defects.

`missing/ci_example` — see its own README.
