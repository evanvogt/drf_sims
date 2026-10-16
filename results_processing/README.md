# Results processing

Figures and summaries built from the `*_metrics.RDS` files each study produces.
Nothing here reads per-run simulation output directly — run the study's
`*_collect.R` and `*_metrics.R` first.

## Layout

| | |
|---|---|
| `thesis_figures/` | the figures that go in the thesis, one script per chapter section |
| `thesis_tables/` | the LaTeX tables that go in the thesis |
| `ICTMC_poster_planning.md` | notes from choosing the ICTMC poster figures (the poster scripts are no longer in the repo) |

Each script reads `../results/<study>/..._metrics.RDS` (for example
`../results/competing_risk/surv_metrics.RDS`) and writes under
`../results/thesis_figures/` or `../results/thesis_tables/`. Run a script from
its own folder: `cd results_processing/thesis_figures && Rscript surv.R`.

## `thesis_figures/`

| script | reads |
|---|---|
| `sample_size.R` | the correlated sample-size studies (the chapter): per outcome, scenarios 1–4 and 5–10 (`*_supp_*`), each at ρ = 0.5 and as the paired ρ = 0.5 − ρ = 0 difference (`*_rhodiff_*`): bias/RMSE, and HTE test rejection rates at `HTE_ALPHA` (`R/figures.R`) |
| `ss_ci.R` | the correlated CI studies, `CI_sf` sweep: per outcome, at ρ = 0, ρ = 0.5 and the paired difference, one figure per interval metric (`{cts,bin}_corr_{rho0,rho05,rhodiff}_ci_{marg,simul,len}.png`), scenarios across, n down, `CI_sf` on x; the inbuilt causal-forest interval is a flat line with its MCSE band. Coverage on a fixed y scale, length free per panel. A study with incomplete metrics is skipped. optimal_sf figures to follow |
| `missing.R` | the continuous missing-data study (main chapter), scenarios 1–4, DR-RandomForest only: bias difference and relative efficiency against the complete-data arm, as separate figures, handling method on x, complete vs incomplete units by colour |
| `miss_cts.R`, `miss_bin.R` | the missing-data studies on the pre-2026-09-29 metrics (all-unit bias/MSE). `missing.R` supersedes `miss_cts.R` for the continuous study; `miss_bin.R` has no equivalent yet |
| `miss_cts_ci.R` | MI confidence intervals |
| `model_eval.R` | the model-evaluation strategies tree (`model_evaluation_strategies/me_strat_metrics.RDS`), the DR / IF proxies with estimated or fixed π: top-1 selection accuracy and regret (`me_top1.png`, `me_regret.png`), with fold regime down, scenario across, n on x, proxy by colour and pipeline by shape; and regret marginalised over each design axis (`me_axis_regret.png`). Regret is drawn against random-pick and best-fixed-choice references. Optional `[metrics.RDS] [figure dir]` arguments for a smoke run |
| `dgm.R` | **no results** — the DGMs' true CATE (appendix), via `../dgm_population.R`, at the primary ρ = 0.5: continuous/binary scenarios 1–4 as τ(x) − ATE (free of n): against X4 by X3 (`dgm_cate_x4.png`), its population distribution with SD(τ) (`dgm_cate_dist.png`), scenario 3 over (X4, X5) by X3 with the covariate density (`dgm_cate_sc3_surface.png`), and against the control mean / risk (`dgm_m0_cate.png`); competing risks, all 7 scenarios, the RMTL CATE of each event by X3 (`dgm_surv_cate.png`). A few minutes locally (numerical integration) |
| `rel_eff_snr.R` | the correlated sample-size studies and the competing-risks study at ρ = 0.5: each run's RMSE(DR) / RMSE(causal forest) against SD(τ) / SD(noise) from the DGM (noise given X and W; for competing risks, of the days lost to the event by the horizon, computed by quadrature). DR-RandomForest and DR-SuperLearner vs Causal forest for continuous/binary; DR-RF on PV, DR-SuperLearner on PV and DR-RSF vs Causal forest on PV for competing risks. `rel_eff_snr.png` (continuous and binary scenarios 1–4, plus competing risks), `rel_eff_snr_supp.png` (scenarios 5–10), and one figure per outcome. Optional `[results root] [figure dir]` arguments for a smoke run |
| `surv.R` | the competing-risks study, the production arm of each estimator family (`SURV_ARM_LABELS`, `R/figures.R`): one figure per metric (CATE bias, RMSE, Pearson) at ρ = 0.5, ρ = 0 and the paired difference (`surv_{rho05,rho0,rhodiff}_{bias,rmse,corr}.png`), censoring on x, model by colour, cause-specific vs subdistribution estimand by shape, event (incl. Combined) down, scenario across |

Shared presentation lives in `R/figures.R`: the display labels
(`MODEL_LABELS`, `METHOD_LABELS`, `MECHANISM_LABELS`), the palette and theme,
`summarise_metrics()` for mean ± Monte Carlo SE, `save_fig()` for the standard
21×15cm size, and the plot shapes both pairs of scripts use
(`point_range_plot`, `distribution_plot`, `make_bm_plots`), so renaming an
estimator is a one-line change.

Each script still owns its paths, its filters and its choice of panels — the
shared module handles presentation, not analysis.

## `thesis_tables/`

Booktabs `longtable`s, written to `../results/thesis_tables/`. For the
sample-size studies, one table for scenarios 1–4 (`*_main.tex`) and one for
5–10 (`*_supp.tex`) per outcome; for the missing-covariate studies, the
appendix's full grid, one table per metric per outcome:

| script | writes | rows × columns |
|---|---|---|
| `ss_tables.R` | `{cts,bin}_corr_{rho05,rhodiff}_ss_{main,supp}.tex` | scenario / n / model × estimation metrics |
| `ss_test_tables.R` | `{cts,bin}_corr_{rho05,rhodiff}_ss_tests_{main,supp}.tex` | scenario / test / model × n, HTE test rejection rates at `HTE_ALPHA` (or their paired ρ difference), with a True CATE reference row |
| `ss_ci_tables.R` | `{cts,bin}_corr_{rho0,rho05,rhodiff}_ci_{marg,simul,len}.tex`, `cts_corr_ci_sf.tex` | correlated CI studies: scenario / n / model × the 10 `CI_sf` values, one interval metric per table (the inbuilt causal-forest interval, which does not depend on `CI_sf`, spans the columns). optimal_sf: scenario / n / ρ (0, 0.5, paired Δ) × ratio picked, plug-in and true-τ coverage, length. A study with incomplete metrics is skipped |
| `miss_tables.R` | `{cts,bin}_miss_<metric>.tex`, `{cts,bin}_miss_tests_<test>.tex` | mechanism / scenario / method × model, one metric (or HTE test, with a True CATE column) per table |
| `surv_tables.R` | `surv_main_{rho05,rho0}.tex`, `surv_<family>_{rho05,rho0}.tex` | competing risks: scenario / censoring / arm × (event × CATE bias, RMSE, Pearson). Main: one production arm per family. Secondary, one per multi-arm family: the raw arms, then paired per-run differences between arms that differ by one factor (fitting or pseudo-values) |
| `me_tables.R` | `me_{top1,rank,regret}.tex`, `me_axis.tex`, `me_fixed.tex`, `me_runs.tex`, `me_if_population.tex`, `me_if_selection.tex` | model evaluation (appendix), the 24 proxies of `thesis_figures/model_eval.R`: scenario / n / score type × fold regime × machine, one measure per table (top-1 accuracy, true rank of the pick, regret), with random-pick and best-fixed-choice rows spanning the columns. Axis: scenario × n / design axis / level × the three measures, each level averaged within run. Fixed: scenario / n / candidate × true PEHE and the three measures if always deployed, the best fixed choice in bold. Runs: the proxy tables' layout, each cell the runs that proxy was scored in, with an all-runs row (full everywhere since the scoring propensity is trimmed, DEVLOG bug S; a short cell means a new failure); every caption's run count is its own table's, read off the data, and the short proxies are listed on the console. IF population: `model_evaluation/me_infl_investigation.R`'s `population.RDS`, scenario / constructed candidate × true PEHE, DR, IF with the true τ, IF with the pseudo-outcome and MSE against μ₁, oracle nuisances, rows by true PEHE and each score's rank in parentheses (skipped with a message if the file is missing). IF selection: CV / XGBoost / π fixed at 0.5, scenario / n × regret under DR, both IF surrogates and a random pick, and the T-learner IF's Spearman correlation with true PEHE and with `mu1_mse`. All portrait. Selection derivation shared with the figures via `model_evaluation/me_selection.R`. Optional `[metrics.RDS] [table dir] [population.RDS]` arguments for a smoke run |
| `dgm_tables.R` | `dgm_cts_bin.tex`, `dgm_surv.tex` | **no results** — the numbers behind `thesis_figures/dgm.R`, from the same draws (`../dgm_population.R`). Continuous/binary: outcome / scenario × SD(τ), 1–99% of τ − ATE, cor(m0, τ), linear R², and per n the ATE with the share of the population whose CATE has the opposite sign. Competing risks: scenario / event × mean, SD, mean by X3, X3's effect with X1, X2 held. Portrait |
| `val_tables.R` | `val_sg_{top,bottom}.tex`, `val_sg_counts.tex`, `val_nom_{tevim,shap}.tex`, `val_sg_se.tex` | interim-analysis validation (appendix to `thesis_figures/validation.R`), HC3 unless stated. Subgroup validated in stage 1's direction: scenario / estimator × the 11 interim points, one table per subgroup, with the runs left after dropping empty stage 2 subgroups (landscape). Counts at 0.25 / 0.5 / 0.75: scenario / estimator / interim × per subgroup runs, significant same / opposite direction, true contrast in stage 1's direction. Stage 1's top covariate: scenario / estimator / interim × X4, X3+X5, X1+X2, X01–X03, X04+X05, one table per importance measure, true modifiers in bold. Classical vs HC3: scenario / subgroup / estimator × landmark × SE kind. Optional `[metrics.RDS] [table dir]` arguments for a smoke run |

Each cell is `mean (MCSE)`, summarised with the figures' `summarise_metrics()`;
the layout helpers are in `R/tables.R`. Run from
`results_processing/thesis_tables/`:

```
Rscript ss_tables.R
Rscript ss_test_tables.R
Rscript ss_ci_tables.R
Rscript miss_tables.R
Rscript surv_tables.R
Rscript me_tables.R
Rscript val_tables.R
Rscript dgm_tables.R
```

The thesis preamble needs `\usepackage{booktabs, longtable, pdflscape, array}`.
Each script's header explains its metric choices (e.g. one bias column,
relative ATE bias but not relative CATE bias).

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
