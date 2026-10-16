# Development log

The history of the project's code and of the re-run campaign. None of it is
needed to read or run the code: it explains why some results were re-run, which
bugs were found along the way, and how the campaign was tracked on the cluster.
The README is the starting point; this file is the record behind it.

Contents: the bug ledger (A–S), the status of each study's re-run (as planned
in September–October 2026), how re-runs are tracked on the cluster, why
resource profiling was removed, the 2026-09-26 scenario renumbering, why `R/`
exists, and the per-study history that used to sit in each README's Status
section.

Other READMEs cite this file for "the bug ledger" and "Resource profiling
(removed)". The ledger, Status and profiling sections are as moved out of the
root README on 2026-10-11.

## Bug ledger

Found during the de-duplication. Each is written up in the relevant folder README.

| | what | where | fixed |
|---|---|---|---|
| A | ran on the **continuous coefficient table** on a logit scale | `confidence_intervals/binary` | yes — re-run |
| B | collect looked for `AUX` where everything else said `MNAR` | `missing/continuous` | yes |
| C | `rel_efficiency` was `NA` everywhere — the reference arm was never collected | `missing/*` | yes |
| D | grid filtered after `expand.grid`, so the array index meant two different things | `missing/*` | yes — structurally |
| E | `.gitignore` `*test*.*` hid the verification harness from git | repo | yes |
| F | `stage_2_sl` discarded the pretested SuperLearner library; the correct branch was dead code | shared | yes — re-run |
| G | `bias` and `ate_bias` had **opposite signs** | all metrics | yes |
| H | propensity trimming only on the SuperLearner path | shared | yes — resolved when `nuisance_rf` moved to whole-sample OOB and picked up `trim_ps` too |
| I | competing risks never validates its SuperLearner library | `competing_risk` | yes — resolved when it adopted the shared crossfitting strategy and picked up `pretest_superlearner`; the split-pseudo T-learner branch is the one remaining unvalidated call site |
| J | stale comments and filenames | various | yes |
| K | `pretest_superlearner()` could return an empty SL library, crashing `nuisance_sl`/`stage_2_sl` | shared | yes |
| L | `run_blp_whole()` had no `tryCatch`, crashed on a constant/degenerate CATE | shared | yes |
| M | `dr_oracle` handed log-odds as outcome predictions — `6b06db3` dropped the `plogis` from the oracle formulas but not `oracle_link = "identity"` | `missing/binary` | yes — no re-run: every finished result predates it |
| N | binary MNAR-Y truth evaluated at U = 0 rather than averaged over U — equal on the identity scale, not the logit one | `missing/binary` | yes — repaired at metrics time from the saved truth; no re-run |
| O | continuous `bW` calibration used `sd = s_err + s2` (ignoring b1, b2 and the heterogeneity variance, and adding SDs) and set `bW` rather than the ATE, so realised power ran 3–100% and scenarios 3, 5, 8 had a positive ATE; b0, b1, b2 also varied by scenario. Now one baseline (0.4, −0.5, 1), and each trial planned for 80% power under homogeneity with the true ATE equal to the planned effect — realised power 61–80% as heterogeneity grows | `R/dgm_scenarios.R` — every continuous study | yes — re-run |
| P | binary `bW` calibration planned at 75% power at the risk plogis(b0), ignoring X1 and X2 (0.401 against a population risk of 0.453), and set `bW` rather than the ATE (a marginal risk difference), so `bW` was the same in every scenario and realised power ran 5–99%; the modifiers in scenarios 2, 5, 6 and 7 also had the opposite sign to the continuous ones, and `b5` was a column no scenario used. Now each trial is planned for 80% power under homogeneity with the true RD equal to the planned effect in every scenario (realised power 79–81%), the signs match `sample_size/continuous/`, and `b5` is gone. *The risk-difference DGM (2026-09-26) has since reversed every binary modifier's sign relative to `sample_size/continuous/`, on purpose — see `sample_size/binary/README.md`* | `R/dgm_scenarios.R` — every binary study | yes — re-run |
| Q | `pretest_superlearner()` dropped a learner on **any warning** (`tryCatch(warning = )` aborts the fit at the first one), so benign glm/gam warnings — and SL.mean's "All algorithms have zero weight" for a propensity near 0.5 — cut folds' libraries to one or two learners; it also kept a learner if even one prediction was non-NA. Now warnings are recorded and muffled, a learner is dropped only on an error or non-finite predictions, and the drops are saved as `sl_dropped`. Fixed alongside the move to per-nuisance libraries (`R/sl_library.R`) | shared, `crossfitting`, `competing_risk` | yes — re-run (SuperLearner arms only) |
| R | the **scoring** nuisance pipelines (XGBoost and AutoML in `me_nuisance.R`) still fit one outcome model on (W, X) and read it at W = 0/1 (an S-learner) after the 9 candidates moved to per-arm outcome models on 2026-10-08, so every φ and every proxy score was S-learner based - an inconsistency with the candidates, not a demonstrated accuracy loss (a local XGBoost check found cor(μ̂₁ − μ̂₀, τ) about the same under both). Now `mu0_DR`/`mu1_DR` are fit per arm in both pipelines and every arm (2026-10-14), for consistency, and `me_strategies.R` recomputes `whole` rather than copying the main tree's | `model_evaluation` | yes — re-run (strategies and split trees) |
| S | the **scoring** pseudo-outcome divided by π̂(1 − π̂) with no trimming, unlike the candidates' (`trim_ps()`). AutoML's propensity was exactly 0 or 1 on some rows of 294/1200 `holdout`, 58 `whole` and 6 `cv_shared` runs (XGBoost never), which made φ NaN or ±Inf and every estimated-π score undefined for the run - summaries on as few as 65% of runs. Where φ was +Inf, the DR score was +Inf for all 9 candidates and the argmin returned the first candidate as a pick. It had been put down to H2O returning NA predictions and excused by a hard-coded list in `me_strategies_verify.R`. Now φ and the IF score use π̂ trimmed to [0.05, 0.95] in every arm and pipeline, `me_selection.R` keeps only finite scores, and verify checks for a finite trimmed φ (2026-10-14) | `model_evaluation` | yes — metrics only: `me_metrics.R` recomputes φ from the saved predictions |

Three more surfaced along the way:

- `missing/binary` was a **half-converted fork**: continuous coefficients, a
  continuous power calculation, and truth on the **log-odds** scale while every
  estimator targets a risk difference
- `binary`'s grid was declared three ways and they disagreed. Submitting
  indices 1–1600 against the analysis script's ten-scenario grid ran runs 1–40
  of all ten scenarios, so the results on disk have 40 replicates per cell,
  not 100 (see `sample_size/binary/bin_config.R` at tag `independent-ss-final`). The study re-runs in full anyway
- `combine_mi()` in `missing/ci_example` read `alpha` as a **free variable** from
  the global environment

## Status

**Retired 2026-10-06:** the independent-covariate sample-size studies
(`continuous`, `binary`, and `confidence_intervals/{continuous, binary,
optimal_sf}`) are not used in the thesis. The correlated studies replace them,
with their ρ = 0 arm as the independent-covariate comparator. Their code is at
git tag `independent-ss-final`, and their READMEs stay for the method notes.
Their results are archived with
`R/archive_old_results.R --label retired_2026-10 --trees continuous binary confidence_intervals`.
The rows below that name them are historical.

**Step 0 — archive the old results.** Nothing under `../results` is what the
current code produces (bug O and the risk-difference DGM, below), and every
`scenario_<k>/` directory there uses the pre-2026-09-26 numbers. Re-running into
those trees would overwrite some results and leave others to be collected under
the wrong scenario, so pack them away first — on the cluster, from the repo
root:

```bash
Rscript R/archive_old_results.R              # dry run: what goes where
cd R/jobscripts; qsub archive_old_results.sh # one .tar per study, then deletes the tree
```

Each study's tree becomes `results/_archive/pre_2026-09-26/<path>.tar`
(uncompressed: the result files are gzipped already), and is deleted only once
its `.tar` lists every file; `tar -xf <that .tar>` from `results/` restores it.
Every study below then re-runs from empty. `competing_risk/` is left alone by this step; its own 2026-10-10 archive is
`Rscript R/archive_old_results.R --label pre_2026-10-10 --trees competing_risk` (see the `competing_risk` row below).

| study | |
|---|---|
| `continuous`, `binary`, `missing/continuous`, `missing/binary`, `missing/ci_example`, `confidence_intervals/continuous`, `confidence_intervals/binary`, `confidence_intervals/optimal_sf` (cts), `confidence_intervals/optimal_sf` (bin), `validation/continuous` | **re-run — crossfitting strategy changed** (see Methods above), on top of any bug-fix re-run already listed below |
| `crossfitting` | its own comparison arms are unchanged by the crossfitting change — but it re-runs for bug O, below |
| `continuous`, `missing/continuous` | also re-run for bug F (`dr_superlearner` only) |
| `binary` | also re-run for bug F (`dr_superlearner` only) |
| `missing/binary` | also re-run — the DGM was wrong three ways |
| `confidence_intervals/binary`, `confidence_intervals/optimal_sf` (bin) | also re-run — the DGM was wrong |
| every continuous study: `continuous`, `missing/continuous`, `missing/ci_example`, `confidence_intervals/continuous`, `confidence_intervals/optimal_sf` (cts), `validation/continuous`, `model_evaluation`, `crossfitting`, `crossfitting/confidence_intervals` | also re-run for bug O — the shared baseline and `bW` changed, which moves every dataset and the level of the true CATE. For `crossfitting`, `crossfitting/confidence_intervals` and `model_evaluation` this is new: any finished results they hold are stale. See `continuous/README.md`'s "Outcome model and `bW` calibration" |
| every binary study: `binary`, `missing/binary`, `confidence_intervals/binary`, `confidence_intervals/optimal_sf` (bin) | also re-run for bug P and the **risk-difference DGM**. Bug P changed the `bW` calibration and the modifier signs in scenarios 2, 5, 6 and 7. The risk-difference DGM then put the treatment effect on the risk-difference scale, with a bounded control risk and rescaled, sign-reversed modifiers. Together they move every dataset's outcome and the true CATE. Submit on the risk-difference code; any bug-P re-run submitted on the logit-scale DGM is superseded too. For `missing/binary` this is new: its finished rows 1–9,900 are superseded, so all 12,600 rows re-run. See `binary/README.md`'s "Outcome model and `bW` calibration" |
| `competing_risk` | **first run under the new strategy** — it has now adopted the crossfitting change (it was the last production study still double-crossfitting) and runs clean end-to-end. Its pseudo-value and SuperLearner frameworks each ship in several arms, because it also crosses a second factor — whole-sample vs crossfit pseudo-values. See its README |
| `model_evaluation` | independent estimator/nuisance code (see its README), whose 9 candidates moved onto the shared single-crossfit strategy. Its 358/360 runs on that strategy predate bug O, so they are archived with the rest and the study re-runs — main, strategies and split trees |

The table above is the September plan (roughly 32,000 array jobs) and is no longer kept up to date. The live count is
`check_all_studies.md`, generated from `R/study_registry.R` (16 studies, about 128,000 array jobs once
`validation/binary` and the correlated studies are included). The 2026-10-10 competing-risks fix (IPW horizon, X4/X5
kept, semi-oracle arms) means every earlier competing-risks result is stale, whatever `check_all_studies.md` says if it
was generated before then. Bug G costs no cluster time: it is computed from the
saved `*_all.RDS` files, so only the metrics scripts and the figures rerun.

### Tracking the rerun

Each study still resubmits its own failures the usual way:
`Rscript <study>/<prefix>_check.R` writes `jobscripts/failed_ids.txt`, then
`qsub jobscripts/<prefix>_rerun.sh` resubmits them. The check script also
rewrites that rerun jobscript - `-J` to match the number of failures, and the
resource request to sit above `<prefix>_1.sh`'s (one more core, 1.2x the
memory, 2x the walltime, capped at 72 hours and never lowered below what the
script already asks for). There is no `-J` to set by hand any more.

Every jobscript writes its PBS output to a `jobscripts/logs*/` directory, and
`*logs*/` is gitignored, so none of them exist after a fresh clone and PBS
rejects the job. Create them all on the cluster before the first submit:

```bash
Rscript make_log_dirs.R            # --dry-run to just list them
```

For a bird's-eye view
across every study at once - how many jobs are expected, found and missing,
and why each study is being rerun - run:

```bash
Rscript check_all.R
```

on the HPC login node (results only exist there - too large to sync back).
It's read-only (never writes `failed_ids.txt`, never calls `qsub`), so it's
safe to run as often as useful during the campaign. It writes
`check_all_studies.csv`/`.md` next to itself; those two files are committed,
so `git push` on the HPC + `git pull` locally is how progress gets checked
from off the cluster. The registry of studies it scans lives in
`R/study_registry.R` - add a row there for any new study.

## Resource profiling (removed)

Profiling with the [`syrup`](https://simonpcouch.github.io/syrup/) package is
**not effective for this study**, and the profiling scripts have been deleted.
Each study used to have a `<prefix>_profile.R` sweep (timing / memory / CPU over
`workers` × `grf_threads`, instrumented with `syrup`), a `jobscripts/<prefix>_profile.sh`
array job, and a `<prefix>_profile_summary.R` that wrote the measured `#PBS -l`
directives into `<prefix>_1.sh`. None of it produced usable numbers:

- on the cluster, `syrup()` starts its sampler with `callr::r_session$new()` and a
  hardcoded 3 s timeout that it gives no way to raise. The profiling subjobs died
  with `Could not start R session, timed out` — most likely because starting an R
  session off the networked `/rds` filesystem takes longer than that
- on cells with `workers > 1` the `future::multisession` workers never showed up
  in the process tree `syrup` samples, so the CPU / memory figures missed them
- locally `syrup` measures nothing useful on Windows

So the `#PBS -l` lines and trailing `Rscript` arguments in the `*_1.sh`
jobscripts are hand-set, not measured. Any that still say *placeholder* stay
that way — adjust them by hand from how runs actually behave (the `_rerun.sh`
resource bump above does part of this automatically). The study READMEs, code
comments and jobscript headers have been updated to match.

The deleted files (the 26 `*_profile*.R`/`.sh` scripts, plus
`crossfitting/cf_diagnose_sampler.R`, `cf_diagnose_multisession.R` and
`jobscripts/cf_diagnose.sh`, which diagnosed the failures above) are in git
history:

```bash
git log --diff-filter=D --name-only -- '*_profile*' '*cf_diagnose*'   # find the deleting commit
git checkout <commit>^ -- crossfitting/cf_profile.R                  # restore a file from its parent
```

## Scenario renumbering (2026-09-26)

The binary and continuous studies share ten scenarios. On 2026-09-26 they were
renumbered so the four the chapters report come first. No data changed: each
scenario generates exactly what it did under its old number. Everything saved
before then uses the old numbers: `../collected_metrics/`, results
directories, figures and the `results_processing/` notebooks.

| now | treatment effect | before 2026-09-26 |
|---|---|---|
| 1 | none (null) | 1 |
| 2 | simple — continuous X4 | 3 |
| 3 | complex — X3 + X4 + X4·X5 | 8 |
| 4 | non-linear — cos(X4) | 9 |
| 5 | binary X3 | 2 |
| 6 | X3 + X4 | 4 |
| 7 | X3·X4 | 5 |
| 8 | X3 + X4 + X3·X4 | 6 |
| 9 | X4·X5 | 7 |
| 10 | exponential in X4 | 10 |

`missing/` had its own numbering before then:

| now | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| `missing/` before 2026-09-26 | 1 | 6 | 4 | 5 | 2 | 3 |

Subsets: `model_evaluation/` ran 1, 4, 6, 8 (`crossfitting/`'s subset) until
2026-10-11, then 1–4; `validation/` ran scenario 2 only until 2026-10-11, then
2–4, and 1–4 from 2026-10-12; `missing/ci_example/` ran 1 and 3–6 until
2026-10-01, then 1–4.

## Why `R/` exists

The repo grew by copy-paste: each new study started as a duplicate of
`sample_size/continuous/`. By the time that stopped, the same CATE estimation
code lived in seven files, the same DGM in four, and the same collect/check
boilerplate in eight; `cts_models.R` and `bin_models.R` differed in two places
out of 438 lines. Consolidating into `R/` removed about 3,000 lines and, more
usefully, the possibility of the copies drifting apart, which is where most of
the bugs in the ledger came from. `crossfitting/` was the model: it already
sourced shared code rather than forking it. The figures went the same way into
`R/figures.R` (`cts_ss.R`/`bin_ss.R` differed in 18 of 280 lines,
`miss_cts.R`/`miss_bin.R` in 64 of 140).

## Study history moved out of the READMEs (2026-10-16)

Each study README used to carry a Status section and dated "since/until"
notes. They are summarised here; the full text is in git at `9f7b43e`, e.g.
`git show 9f7b43e:model_evaluation/README.md`.

**`model_evaluation/`.** Ported from a single-commit prototype that never
completed a run and used the external `benchtm` package for data; the port
fixed nonexistent `here("src", ...)` paths, an undefined variable, a `future`
plan leak, a `NUL` typo, missing `on.exit()` parens, a dead `fix_automl()`
check, a duplicate `collate_predictions()`, a missing `n_cores` argument and
`rf2`/`rf3`'s `mtry` (calibrated to `benchtm`'s wider covariate set). The
candidates moved from double to single crossfitting after the port. 2026-10-08:
correlated covariates (ρ = 0.5), 30 → 100 runs, per-arm (T-learner) outcome
models in all 9 candidates (the S-learner made `net*`'s `Y1.hat − Y0.hat`
constant), mean + glm propensity for `SL*`, `SL1` = the production
`dr_superlearner`, `net2` ridge → interaction lasso, and the `cv_indep` arm
(a second, independent fold draw) removed. 2026-10-11: scenarios 1, 4, 6, 8 →
1–4, listed as `c(1, 4, 2, 3)` so 1 and 4 kept their array indices.
2026-10-14: per-arm scoring nuisances (bug R) and trimmed scoring propensity
(bug S); untrimmed, AutoML's propensity was exactly 0 or 1 on some rows of
294/1200 `holdout`, 58 `whole` and 6 `cv_shared` runs. 2026-10-15: `inflt`,
`inflt05` and `mu1_mse` added (metrics only). The pre-bug-O tree had 358/360
main runs, 69 of 358 strategies runs with NA automl `holdout` phi.

**`competing_risk/`.** Double crossfitting replaced by the shared strategy;
outcome models S-learner → T-learner on 2026-09-27; trimmed propensities and
pretested SuperLearner libraries came with it. 2026-10-01: event 1 baseline
scale 15 → 25 and log-HR-matched effects (the old values ran the scenarios
together on the RMTL scale); correlated covariates, ρ ∈ {0, 0.5}; the RSF
DR-learner arms. 2026-10-06: ρ = 0.5 primary, ρ = 0 sensitivity (the reverse
before). 2026-10-10: the copula's X4/X5 kept as noise covariates (previously
drawn and dropped); `get_ipw()` stopped dropping units censored or with the
competing event after the horizon (`ipw` was badly biased, event 2 worst);
semi-oracle DR arms. On the pre-2026-10-01 DGM a full rerun left 19
deterministic failures: 15 from whole-sample pseudo-values NA when nobody was
observed past the horizon, 4 from `pretest_superlearner()`'s SL.mean fallback
on a constant treated-arm RMTL2 cell; an intermediate 398 were probably
walltime (the jobscript had been cut to 1 h). The PBS logs for that rerun were
deleted before anyone read them; the jobscripts have used `#PBS -j oe` since.

**`validation/`.** Rebuilt from a pre-restructure fork (standalone DGM and
estimator copies, hand-written grid arithmetic in `validation.sh`, hardcoded
cluster paths, a broken collect/metrics handoff) onto `R/`; the causal-forest
nuisance changed to `Y.hat.cf`. `bottom_pval` read `pvals_bottom[1]`, the
intercept, so it was never a subgroup test; `interaction_pval()` now indexes
`W:v` by name. Interim grid 3 → 11 points; TreeSHAP and `top_var_tests` added.
2026-10-08: correlated covariates, one trial split rather than two draws,
`p_cts_adj`, the DR SuperLearner, the code moved to `val_common.R` for the new
binary arm, HC3 interaction tests. 2026-10-09: ρ = 0 briefly in the grid.
2026-10-11: scenarios 2–4, classical and HC3 SEs both saved, `p_dd`/`p_proj`,
`median_split()` fixed for a binary covariate (it was NA whenever X3 was
nominated), per-scenario reports. 2026-10-12: null scenario 1 (4,400 jobs,
grid reordered), interaction estimates and `ref1`/`true2` contrasts.

**`missing/`.** Bugs B, C, D, M and N (ledger above) were found here.
`missing/binary` was a half-converted fork of the continuous study. 2026-09-27:
the MI arm's per-imputation tests saved (earlier MI runs had none, and no
nuisances to recover them). 2026-09-28: correlated covariates with X01–X03 as
auxiliaries; mechanisms MAR / MNAR / MNAR-Y (and `ci_example`'s `AUX`/`AUX-Y`)
replaced by MAR / MNAR-Y0 / MNAR-tau — the old "MNAR" had `U` in neither X nor
Y, so was MCAR in effect, and the old "MNAR-Y" is MNAR-tau; `ci_example`'s
`mib`/`hybrid` critical value `1 - alpha/2` → `1 - alpha` (97.5% bands
before). 2026-09-29: `U` drawn under MAR (paired mechanisms), missingness mask
saved, split `_cu`/`_iu` scoring, `missforest`/`multiple_imputation` within arm
with Y; binary `bU` 0.08 → 0.12, `sU` 1 → 2 and its own `RD_SCALE_MISS`
(scenario 4 0.204 → 0.179). 2026-09-30: `IPW` model `+ W * Y` with stabilised
weights, `regression` within arm with Y, MI keeps `tau_mi` and the variance
components, `blp_se_hc3`/`blp_p_os` in `mi_tests`. 2026-10-11: `ci_example`'s
`hybrid` margin was `sqrt(lambda_hat * S_star)`, about half the correct width.
`dr_random_forest` once carried no HTE tests under `profile = "missing"`
(copy-paste drift); the archived results were back-filled by a one-off patch,
removed after commit `e7b1d59`.

**`crossfitting/`.** Its conclusions were first drawn on the pre-bug-O DGM.
2026-09-27: arms rebuilt around per-arm outcome models; the S-vs-T arms
`scf_oob_t`, `oob_oob_s` and `oob_oob_manual` removed. The CI pilot replaced a
`run_crossfit_structured_arms()` orchestrator that trimmed nuisance fits and so
drew different `dcf`/`cf_dcf` values from the production study.
