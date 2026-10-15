# Correlated binary: signal-to-noise check

Does the DR-learner close the gap on the causal forest when there is more
signal? `results_processing/thesis_figures/rel_eff_snr.R` plots the binary
DR-learners against the causal forest by SD(τ)/SD(noise). In the main study
(`../`, ρ = 0.5) that ratio is at most 0.2 in scenarios 1–4. The binary noise
is Bernoulli, about 0.46 whatever the DGM, so SNR can only rise through
SD(τ). This check moves away from the main study in two directions, one at a
time, as `competing_risk/snr/` does: a larger HTE (higher SNR) and a larger n
(more data at the same SNR).

| | |
|---|---|
| scenarios | 1–4 (null, simple, complex, non-linear), from `../` |
| set | `binary_snr`: `binary_corr_0.5` plus its own risk thresholds (below) |
| ρ | 0.5 (the main study's primary) |
| bW | calibrated at n = 500 in every cell, so the ATE is −0.118 throughout |
| runs | 100 per cell |
| array | **2,500 jobs**: 1–1300 `bin_snr_1.sh` (n = 500), 1301–2500 `bin_snr_2.sh` (n ≥ 1000) |
| results | `../results/correlated/binary_snr/scenario_<k>/<n>/hte_<mult>/` |

## Design

The cells are one factor at a time from the baseline (n = 500, ×1), not a
factorial grid. Run varies fastest, so each cell is a block of 100 rows:

| rows | n | scenarios | `hte_mult` |
|---|---|---|---|
| 1–100 | 500 | 1 | 1 |
| 101–1300 | 500 | 2, 3, 4 (in that order) | 1, 1.5, 2, 3 within each scenario |
| 1301–2500 | 1000, 2000, 4000 (in that order) | 1–4 within each n | 1 |

So index 101 is scenario 2 ×1 run 1, 1201 is scenario 4 ×3 run 1, and 2101 is
n = 4000 scenario 1 run 1. Scenario 1 has no HTE, so it is on the n axis only.

Population SNR is `rel_eff_snr.R`'s definition, SD(τ) / sqrt(E[p(1 − p)]),
computed on the clamped truth. These figures are from `R/bin_verify_hte.R`
check 8 (10^5 draws); `bin_snr_metrics.R` recomputes them on 10^6 draws:

| scenario | ×1 | ×1.5 | ×2 | ×3 | clamped at ×3 | ATE at ×3 |
|---|---|---|---|---|---|---|
| 2 (tanh X4) | 0.111 | 0.167 | 0.224 | 0.340 | 3.2% | −0.118 |
| 3 (X3 + X4 + X4·X5) | 0.120 | 0.181 | 0.243 | 0.371 | 0% | −0.118 |
| 4 (cos X4) | 0.196 | 0.298 | 0.404 | 0.606 | 6.7% | −0.123 |

SNR is a population quantity and bW is fixed, so the n cells sit at their
scenario's ×1 SNR. Read them against n.

- **`hte_mult` scales the effect modifiers** b3, b4, b34 and b45, through
  `generate_scenario_data(hte_mult = )` and `get_oracle_info(hte_mult = )`.
  bW is then recalibrated, so the ATE is held. That is why scenario 3's and 4's
  bW move with the multiplier: E[g] is not zero there. The default of 1 leaves
  every other study bit-identical (`tests/run_tests.R`).
- **The risk thresholds diverge from the main study's.** Under the main rule
  (the worst-case treated risk over the covariate support inside
  [0.01, 0.99]), the multipliers would stop at about ×2.6, ×2.5 and ×1.25 for
  scenarios 2, 3 and 4. This set clamps the treated risk into [0.01, 0.99]
  instead, and the truth is the **clamped** risk difference (`truth_at()`, and
  the oracle formula through `get_oracle_info()`). The control risk lies in
  [0.34, 0.70], so only the treated risk is ever clamped. Where it is, τ(x) is
  no longer exactly the scenario's `te_expr`, and X1 and X2 modify it a
  little. At ×1 nothing is clamped. `calibrate_bW()` ignores the clamp, so the
  ATE drifts slightly in the clamped cells (−0.123 against −0.118 at
  scenario 4 ×3).
- **bW is calibrated at n = 500 in every cell** (`calib_n`). The main study
  recalibrates per n, giving a smaller ATE at larger n. Here the n cells share
  one DGM, so only the amount of data changes.
- **Seeding** is `setup_rng_stream(run)`, as in the main study. The n = 500
  cells of a run share W and X across scenarios and multipliers, so the
  `hte_mult` comparison is paired. Y is coupled only through rbinom's
  uniforms. The n cells are not paired.
- **The baseline equals the main study.** The SuperLearner arm runs last in
  `cate_methods()`, so leaving it out does not move the RNG stream of the arms
  before it. Index 101 reproduces `../bin_corr_analysis.R 8010` (scenario 2,
  n = 500, run 1, ρ = 0.5) bit for bit: the data, the truth,
  `nuisances_rf` and all four arms' τ. This was checked locally on
  2026-10-15.

## Arms

`cate_methods(models = SNR_MODELS, family = binomial(), profile = "base")`,
with no SuperLearner library:

| arm | |
|---|---|
| `causal_forest` | the comparator |
| `dr_random_forest` | DR-learner, grf per-arm outcome forests and propensity |
| `dr_semi_oracle` | DR-learner with the true propensity 0.5 |
| `dr_oracle` | DR-learner with the true (clamped) outcome model and propensity: the stage-2 floor |
| `t_random_forest` | derived at metrics time from `dr_random_forest`'s stage 1 (`add_t_learners()`) |

Unlike competing_risk/snr's semi-oracle twins, `dr_semi_oracle` refits its own
outcome forests (`t_learner_rf()`) rather than reusing `nuisance_rf`'s. Its
gap to `dr_random_forest` is therefore the propensity plus forest Monte Carlo
noise. The "base" profile keeps the BLP and independence tests, so test power
against SNR comes as a by-product.

## Running it

On the cluster, from `sample_size/correlated/binary/snr/jobscripts/` (run
`Rscript make_log_dirs.R` from the repo root first):

```bash
qsub bin_snr_1.sh                 # 1-1300, n = 500
qsub bin_snr_2.sh                 # 1301-2500, n = 1000 / 2000 / 4000
Rscript ../bin_snr_check.R        # writes failed_ids.txt, points bin_snr_rerun.sh at it
qsub bin_snr_rerun.sh             # only if the check found failures
qsub bin_snr_collect.sh
qsub bin_snr_metrics.sh
```

Local smoke test, from `sample_size/correlated/binary/snr/`:
`Rscript bin_snr_analysis.R 101` (scenario 2 ×1, the main DGM) or
`Rscript bin_snr_analysis.R 1201` (scenario 4 ×3, the most clamping). The
arguments are `<i> [grf_threads]`. There is no `workers` argument, because
only the SuperLearner arm uses a future plan.

`bin_snr_metrics.R` writes three files:

- `bin_snr_population.RDS`: per (scenario, `hte_mult`), bW, the ATE, SD(τ),
  SD(noise), SNR and the share clamped.
- `bin_snr_metrics.RDS`: one row per scenario × n × `hte_mult` × run × model,
  with `cate_metrics()` and `hte_test_metrics()`.
- `bin_snr_ratios.RDS`: each run's RMSE of every other arm divided by the same
  run's `causal_forest` RMSE, summarised per cell (median, IQR, geometric mean
  with the MCSE of the mean log ratio, and the share of runs where the arm did
  better), with the cell's SNR.

## Resources

Measured locally on R 4.5.3 with 2 grf threads. The cluster runs grf 2.4.0
against 2.5.0 here, so treat these as rough guides:

| n | wall | peak memory |
|---|---|---|
| 500 | about 10 s | well under 1gb |
| 4000 | about 1.3 min | under 1gb |

There is no rfsrc here, so memory stays small. The jobscripts ask for 2gb
(n = 500) and 4gb (n ≥ 1000), with 15 and 30 minutes, as headroom.
