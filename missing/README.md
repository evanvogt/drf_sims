# Missing covariates

What happens to CATE estimation when covariates are missing, and which handling
method costs least?

Three studies:

| folder | question |
|---|---|
| `continuous/` | the main design — continuous outcome, 8 handling methods + a `complete_data` reference × 3 mechanisms (MAR, MNAR-Y0, MNAR-τ) |
| `binary/` | the same design with a binary outcome, the treatment effect on the risk-difference scale |
| `ci_example/` | how to build a *confidence interval* after multiple imputation |

The full design — the correlated covariates, the mechanisms and their pairing,
the amputation, each handling method, the split scoring on complete and
incomplete units, and the open decision on pooling the MI arm's HTE tests — is
in [`ADEMP.md`](ADEMP.md). Each study folder's README covers its own run
commands and files.

Scenario `k` here is scenario `k` in `sample_size/`. Results saved before
2026-09-26 use older scenario numbers; the map is in `docs/DEVLOG.md`.

## Running

The continuous and binary grids are **12,600 jobs** each: rows 1–9,900 are
scenarios 1, 3, 4, 5 (scenario 1 has no MNAR-τ, so 900 rows are dropped), and
rows 9,901–12,600 are scenario 2, appended later as its own block so the
earlier rows keep their meaning.

```bash
qsub missing/continuous/jobscripts/cts_miss_1.sh       # rows 1-9,900
qsub missing/continuous/jobscripts/cts_miss_extra.sh   # rows 9,901-12,600, scenario 2
qsub missing/binary/jobscripts/bin_miss_1.sh
qsub missing/binary/jobscripts/bin_miss_extra.sh
```

Then check, collect and metrics as usual (see each study's README). Every
model carries its HTE tests natively (`PROFILES$missing`), and the MI arm
saves its per-imputation tests.

`Rscript missing/miss_dgm_checks.R` (from the repo root) checks the DGM's
claims at n = 20,000 with no simulation output.

## Results files

`miss_results.R` / `miss_results.qmd` (one copy for both outcomes: an
`outcome` argument / Quarto param) show **every** metric a study's
`*_metrics.R` computes — bias, ATE bias, MSE, RMSE, MAE, relative
efficiency, correlation, Spearman, sign accuracy, the three HTE-test p-values,
and the `n_na` diagnostic. Figures go to `results/all_figures/missing/<study>/`.

These are the diagnostic view. `results_processing/thesis_figures/miss_*.R` is
the chapter view: the two or three panels that get printed, written to
`results/thesis_figures/`. The two never overwrite each other.

Labels, palette, the mean ± MCSE summary and the point-range panel are shared
via `R/figures.R` — `MISS_SCENARIO_LABELS`, `METHOD_LABELS`,
`MECHANISM_LABELS`, `STRATEGY_LABELS` and `point_range_plot()` in particular.
Rename an estimator or a handling method there and every figure follows.
