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
| scenario | 7 (HTE on both events), from `../surv_dgm.R` |
| ρ | 0.5 (the main study's primary) |
| censoring | TRUE |
| runs | 100 per cell |
| array | **700 jobs**: 1–400 `snr_1.sh`, 401–600 `snr_2.sh`, 601–700 `snr_3.sh` |
| horizon | 28 |
| results | `../results/competing_risk_snr/<n>/hte_<mult>/` |

## Design

The cells are one factor at a time from the baseline, not a factorial grid:

| rows | n | `hte_mult` | E1 HR, X3 = 0 / 1 | E2 HR, X3 = 0 / 1 | SNR Event 1 / Event 2 |
|---|---|---|---|---|---|
| 1–100 | 500 | 1 (baseline) | 1.42 / 2.86 | 0.70 / 0.35 | 0.29 / 0.13 |
| 101–200 | 500 | 1.5 | 1.42 / 4.06 | 0.70 / 0.25 | 0.43 / 0.17 |
| 201–300 | 500 | 2 | 1.42 / 5.75 | 0.70 / 0.17 | 0.57 / 0.20 |
| 301–400 | 500 | 3 | 1.42 / 11.6 | 0.70 / 0.09 | 0.79 / 0.23 |
| 401–500 | 1000 | 1 | as baseline | | 0.29 / 0.13 |
| 501–600 | 2000 | 1 | | | |
| 601–700 | 4000 | 1 | | | |

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
- **Seeding** is `setup_rng_stream(run)`, as in the main study. The four
  n = 500 cells of a run share W, X and U, so the `hte_mult` comparison is
  paired. The n cells are not. The baseline is rerun here rather than reused
  from the main array, because running a subset of arms moves the RNG stream
  of every arm fitted after the ones left out.

## Arms

The production crossfitting arms (README "Crossfitting"), via
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

## Running it

On the cluster, from `competing_risk/snr/jobscripts/` (run
`Rscript make_log_dirs.R` from the repo root first):

```bash
qsub snr_1.sh                 # 1-400, n = 500
qsub snr_2.sh                 # 401-600, n = 1000 / 2000
qsub -J 601-602 snr_3.sh      # n = 4000: two first, see Resources
qsub -J 603-700 snr_3.sh      # then the rest (command-line -J overrides the script's)
Rscript ../snr_check.R        # writes failed_ids.txt, points snr_rerun.sh at it
qsub snr_rerun.sh             # only if the check found failures
qsub snr_collect.sh
qsub snr_metrics.sh
```

Local smoke test, from `competing_risk/snr/`: `Rscript snr_analysis.R 1`
(the baseline, run 1, about 30 s) or `Rscript snr_analysis.R 301` (×3, run 1).

`snr_metrics.R` writes `snr_metrics.RDS` (one row per n × `hte_mult` × run ×
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
n = 4000 cell needs around 25–30gb. `snr_3.sh` asks for 64gb and 2h as
headroom, since that number is not measured. Submit one or two indices first
and read the `.o` file before submitting the rest. `snr_rerun.sh` is sized like
`snr_3.sh`, because a rerun list can mix every n.
