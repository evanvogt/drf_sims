# Competing risks: signal-to-noise check

Does the DR-learner close the gap on the causal forest when there is more
signal? `results_processing/thesis_figures/rel_eff_snr.R` plots the
competing-risk DR-learners against the pseudo-value causal forest by
SD(τ)/SD(noise). In the main study every scenario × event cell is at or below
0.29, and most are below 0.05. This check moves away from that in two
directions, one at a time: a larger X3 effect (higher SNR) and a larger n
(more data at the same SNR).

| | |
|---|---|
| scenario | 1–7 from `../surv_dgm.R`; 7 (HTE on both events) is the baseline the check started from |
| ρ | 0.5 (the main study's primary) |
| censoring | TRUE |
| runs | 100 per cell |
| array | **4,300 jobs**: scenario 7 is 1–700 (1–400 `snr_1.sh`, 401–600 `snr_2.sh`, 601–700 `snr_3.sh`); scenarios 1–6 are 701–4300 (701–2500 `snr_4.sh`, 2501–3700 `snr_5.sh`, 3701–4300 `snr_6.sh`) |
| horizon | 28 |
| results | `../results/competing_risk_snr/scenario_<k>/<n>/hte_<mult>/` |

## Design

The cells are one factor at a time from the baseline (n = 500, `hte_mult` 1),
not a factorial grid, in each scenario. Scenarios 1 and 2 have no X3 effect
modifier (`b3_1 = b3_2 = 0`), so `hte_mult` changes nothing there and they get
the n axis only: 4 cells. Scenarios 3–7 get all 7 cells below. Scenario 7 is
shown; the others are the same cells with their own treatment effects (SNR in
the population table, `snr_population.R`).

| scenario 7 rows | n | `hte_mult` | E1 HR, X3 = 0 / 1 | E2 HR, X3 = 0 / 1 | SNR Event 1 / Event 2 |
|---|---|---|---|---|---|
| 1–100 | 500 | 1 (baseline) | 1.42 / 2.86 | 0.70 / 0.35 | 0.29 / 0.13 |
| 101–200 | 500 | 1.5 | 1.42 / 4.06 | 0.70 / 0.25 | 0.43 / 0.17 |
| 201–300 | 500 | 2 | 1.42 / 5.75 | 0.70 / 0.17 | 0.57 / 0.20 |
| 301–400 | 500 | 3 | 1.42 / 11.6 | 0.70 / 0.09 | 0.79 / 0.23 |
| 401–500 | 1000 | 1 | as baseline | | 0.29 / 0.13 |
| 501–600 | 2000 | 1 | | | |
| 601–700 | 4000 | 1 | | | |

The grid, by block (run varies fastest, 100 runs per cell):

| rows | scenarios | n | `hte_mult` | jobscript |
|---|---|---|---|---|
| 1–400 | 7 | 500 | 1, 1.5, 2, 3 | `snr_1.sh` |
| 401–600 | 7 | 1000, 2000 | 1 | `snr_2.sh` |
| 601–700 | 7 | 4000 | 1 | `snr_3.sh` |
| 701–900 | 1, 2 | 500 | 1 | `snr_4.sh` |
| 901–2500 | 3, 4, 5, 6 | 500 | 1, 1.5, 2, 3 | `snr_4.sh` |
| 2501–3700 | 1–6 | 1000, 2000 | 1 | `snr_5.sh` |
| 3701–4300 | 1–6 | 4000 | 1 | `snr_6.sh` |

Scenario 7 kept rows 1–700 when the other scenarios were added, so its row
numbers did not move. Its results folder did gain a `scenario_7/` level.

SNR is `rel_eff_snr.R`'s definition, population values at ρ = 0.5. SNR is a
population quantity, so the n cells all sit at one SNR. Read them against n.

- **`hte_mult` scales only the X3 effect modifiers** (`b3_1`, `b3_2`), through
  `generate_surv_data(hte_mult = )`. Its default of 1 leaves the main study
  bit-identical (checked against HEAD). The ATE terms are left alone on
  purpose. Scaling them as well gave less SNR, made Event 2's SNR *fall*
  (RMTL2 saturates once nearly every treated unit has event 1), and left about
  0.5% of the treated arm event-free at 28 at ×3. That last part reopens route
  C (README "Known issues": NA whole-sample pseudo-values). With `b3` alone, ×3
  still leaves about 2% of the treated arm event-free at 28. That is about 4
  treated units per run past 28 with censoring on, so expect roughly 1 run in
  100 to fail that way in the ×3 cell.
- **Seeding** is `setup_rng_stream(run)`, as in the main study. The n = 500
  cells of a run within one scenario share W, X and U, so the `hte_mult`
  comparison is paired. The n cells are not. The baseline is rerun here rather than reused
  from the main array, because running a subset of arms moves the RNG stream
  of every arm fitted after the ones left out.

## Arms

The production crossfitting arms (`../ADEMP.md`, "Pseudo-values"), via
`all_cate_surv_models(models = SNR_MODELS, estimands = c("RMTL1", "RMTL2"))`:

| arm | |
|---|---|
| `pseudo_cf_whole_oob` | causal forest on pseudo-values, the comparator |
| `pseudo_dr_whole_oob` | DR-learner on pseudo-values, grf outcome forests |
| `rsf_dr_oob` | DR-learner, competing-risks RSF outcome model |
| `pseudo_dr_whole_oob_semi`, `rsf_dr_oob_semi` | the same two with the true propensity 0.5 in the pseudo-outcome |

There are no SuperLearner, scf/cvps, ipw or csf arms. The semi-oracle twins
reuse their parent's outcome models, so each pair differs only in the
propensity. They show how much of any DR gap comes from estimating e(x). RMSTc
is not fitted, because `rel_eff_snr.R` has no SNR for it.

### Oracle add-on

`pseudo_dr_oracle` is the DR-learner with the true μ₀(x), μ₁(x) (the per-unit
`RMTL*_0` / `RMTL*_1` in `truth`) and e(x) = 0.5, and the same whole-sample OOB
stage 2. One arm covers both DR-learners: they share the pseudo-values in the
correction term and differ only in μ̂. It is still a pseudo-value DR-learner,
so it removes nuisance error but keeps the pseudo-value noise and the stage-2
error. Production vs semi-oracle is the cost of estimating e(x), semi-oracle vs
oracle the cost of estimating μ, and oracle vs CF what is left with both known.

It is fitted after the fact from the saved runs, not in `snr_analysis.R`: each
`res_sim_<run>.RDS` holds `data` and `truth`, and `pseudo_all()` draws no RNG,
so nothing is rerun and there is no RSF. `snr_oracle_config.R` is its own study
on the same grid (row i is row i of `snr_config.R`), with results in
`../results/competing_risk_snr_oracle/` and only the arm's τ̂ in each file. The
stage-2 forests are seeded with `setup_rng_stream(run)`. `snr_metrics.R` merges
`snr_oracle_all.RDS` into each run when it exists, so the oracle appears in
`snr_metrics.RDS` and `snr_ratios.RDS` like any other arm.
`snr_results.qmd` does not plot it yet (`ARMS` leaves it out).

Checked on 4 × n = 5,000 draws of scenario 7: the oracle pseudo-outcome
regressed on the true τ has slope 1.08 (SE 0.05) for Event 1 and 1.09 (0.11)
for Event 2, and the mean of po − τ is within one SE of 0. Its SD is 7× SD(τ)
for Event 1 and 16× for Event 2, so the pseudo-value noise dominates the
nuisance error at these SNRs.

## Running it

On the cluster, from `competing_risk/snr/jobscripts/` (run
`Rscript make_log_dirs.R` from the repo root first):

```bash
qsub snr_1.sh                 # 1-400, scenario 7, n = 500
qsub snr_2.sh                 # 401-600, scenario 7, n = 1000 / 2000
qsub snr_4.sh                 # 701-2500, scenarios 1-6, n = 500
qsub snr_5.sh                 # 2501-3700, scenarios 1-6, n = 1000 / 2000
qsub -J 601-602 snr_3.sh      # scenario 7, n = 4000: two first, see Resources
qsub -J 603-700 snr_3.sh      # then the rest (command-line -J overrides the script's)
qsub -J 3701-3702 snr_6.sh    # scenarios 1-6, n = 4000: two first
qsub -J 3703-4300 snr_6.sh    # then the rest
Rscript ../../../study.R check competing_risk/snr/snr_config.R        # writes failed_ids.txt, points snr_rerun.sh at it
qsub snr_rerun.sh             # only if the check found failures
qsub snr_oracle_1.sh          # oracle add-on, once the main runs are all in
Rscript ../../../study.R check competing_risk/snr/snr_oracle_config.R # writes failed_ids_oracle.txt
qsub snr_oracle_rerun.sh      # only if that check found failures
qsub snr_collect.sh
qsub snr_oracle_collect.sh
qsub snr_metrics.sh           # after both collects
```

`snr_oracle_1.sh` steps its array by 100 (`-J 1-4300:100`), so each of its 43
subjobs fits one cell: `Rscript snr_oracle_analysis.R <first> 100 2`. A row
whose main result file is missing is skipped and shows up in the oracle check.

Local smoke test, from `competing_risk/snr/`: `Rscript snr_analysis.R 1`
(scenario 7 baseline, run 1, about 30 s), `Rscript snr_analysis.R 301` (×3,
run 1) or `Rscript snr_analysis.R 701` (scenario 1, run 1). Each writes a
result under `../results/competing_risk_snr/`. After `snr_analysis.R 1`,
`Rscript snr_oracle_analysis.R 1` fits the oracle on it (about 2 s).

`snr_population.R` writes `snr_population.RDS`, the population SNR of every
scenario × `hte_mult` cell (no cluster, about 15 s).

`snr_metrics.R` writes `snr_metrics.RDS` (one row per scenario × n × `hte_mult` × run ×
arm × event: `cate_metrics()` plus Pearson, Spearman and C) and
`snr_ratios.RDS`. The ratios file holds each run's RMSE of every DR arm divided
by the same run's `pseudo_cf_whole_oob` RMSE, summarised per cell (median, IQR,
geometric mean with the MCSE of the mean log ratio, and the share of runs
where DR did better).

## Resources

Measured locally on R 4.5.3, 2 workers. The cluster runs different package
versions (rfsrc 3.2.2 vs 3.9.0), so treat these as rough guides:

| n | wall | peak memory |
|---|---|---|
| 500 | 0.5 min | 1.3gb |
| 1000 | 1.2 min | 2.7gb |
| 2000 | 3.9 min | 8.0gb |
| 4000 | — | > 9.4gb: rfsrc failed to allocate memory after 8.6 min on a 15gb laptop |

Nearly all of the memory and over half the time at n = 2000 goes on the
competing-risks RSF (`rsf_fit`, `ntime = 0`). It keeps every unique event time
in every terminal node, so its memory grows roughly as n². Extrapolating, the
n = 4000 cell needs around 25–30gb. `snr_3.sh` and `snr_6.sh` ask for 32gb and
2h, since that number is not measured. Submit one or two indices first and read
the `.o` file before submitting the rest. `snr_rerun.sh` is sized like them,
because a rerun list can mix every n.

The timings and memory above are scenario 7's. The other scenarios have not
been timed. `snr_4.sh` and `snr_6.sh` copy `snr_1.sh` and `snr_3.sh`'s
resources, and `snr_5.sh` asks for 12gb rather than `snr_2.sh`'s 8gb, since
8gb is the n = 2000 peak.
