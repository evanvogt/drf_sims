# Results processing

Figures and summaries built from the `*_metrics.RDS` files each study produces.
Nothing here reads per-run simulation output directly — run the study's
`*_collect.R` and `*_metrics.R` first.

## Layout

| | |
|---|---|
| `thesis_figures/` | the figures that go in the thesis, one script per chapter section |
| `LSR_figures.R` | figures for the LSR presentation |
| `results_10_25.Rmd`, `results_update.Rmd`, `bin_new_metrics.Rmd` | working notebooks, snapshots of earlier analyses |

## `thesis_figures/`

| script | reads |
|---|---|
| `cts_ss.R`, `bin_ss.R` | the sample-size studies |
| `miss_cts.R`, `miss_bin.R` | the missing-data studies |
| `cts_ci.R` | continuous confidence intervals |
| `miss_cts_ci.R` | MI confidence intervals |
| `surv.R` | competing risks |

Shared presentation lives in `R/figures.R`: the display labels
(`MODEL_LABELS`, `METHOD_LABELS`, `MECHANISM_LABELS`), the palette and theme,
`summarise_metrics()` for mean ± Monte Carlo SE, `save_fig()` for the standard
21×15cm size, and the plot shapes both pairs of scripts use
(`point_range_plot`, `distribution_plot`, `make_bm_plots`).

That module exists because the scripts had drifted into near-copies:
`cts_ss.R` and `bin_ss.R` differed in 18 of 280 lines once the outcome prefix was
normalised, `miss_cts.R` and `miss_bin.R` in 64 of 140, and `make_bm_plots()` was
defined **identically** in both of the latter. Renaming an estimator for the
write-up used to mean editing six files.

Each script still owns its paths, its filters and its choice of panels — the
shared module handles presentation, not analysis.

## Regenerate after a metrics change

These scripts read `bias` directly with no sign compensation. Bug G flipped the
convention to `est - true`, so **every figure showing bias needs regenerating**,
though none of the scripts needed editing.

`MECHANISM_LABELS` maps both `MNAR`/`MNAR-Y` and the older `AUX`/`AUX-Y`
spellings to the same display label, so figures built from results generated
before the rename still work.

## Gotchas

`ggridges`, `patchwork`, `scales` and `purrr` are loaded by the individual
scripts that use them, not by `R/figures.R` — a shared module should not force
packages on callers that do not need them. If a script fails on a missing
plotting package, that is why.

The `.Rmd` notebooks are historical. They read result files that may predate the
current metric definitions, and are not maintained.
