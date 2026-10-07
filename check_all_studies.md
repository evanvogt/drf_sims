# Rerun campaign status

Last updated: 2026-10-06 18:13 BST

| study_name | category | expected_jobs | found_jobs | missing_jobs | pct_complete | status | reason |
|---|---|---|---|---|---|---|---|
| validation/continuous | crossfit_rerun |  1100 |     0 | 1100 |   0.0 | not_started | crossfitting strategy change; plus T-learner (per-arm) DR outcome models |
| model_evaluation | dgm_rerun |   360 |     0 |  360 |   0.0 | not_started | its 358/360 runs (the first under single crossfitting, me_models.R) predate bug O, which changed the continuous DGM, and the 2026-09-26 renumbering (1/4/6/9 -> 1/4/6/8); archived with the strategies and split trees, so the count restarts from zero |
| missing/ci_example | crossfit_rerun |   400 |   143 |  257 |  35.8 | in_progress | crossfitting strategy change; plus T-learner (per-arm) DR outcome models; plus the 2026-09-28 correlated covariates of continuous_missing (missing/ADEMP.md) |
| correlated/confidence_intervals/optimal_sf (bin) | first_run |  1600 |   274 | 1326 |  17.1 | in_progress | new 2026-10-04: confidence_intervals/optimal_sf (bin) on correlated/'s scenarios 1-4, rho = 0 and 0.5, n 500/1000, 100 runs (sample_size/correlated/confidence_intervals/README.md) |
| correlated/confidence_intervals/optimal_sf (cts) | first_run |  1600 |   259 | 1341 |  16.2 | in_progress | new 2026-10-04: confidence_intervals/optimal_sf (cts) on correlated/'s scenarios 1-4, rho = 0 and 0.5, n 500/1000, 100 runs (sample_size/correlated/confidence_intervals/README.md) |
| competing_risk | crossfit_rerun | 14000 | 14000 |    0 | 100.0 | complete | crossfitting strategy change - the last production study still double-crossfitting; now runs clean end-to-end, so this is its first run under the new strategy; plus bug Q and the per-nuisance SuperLearner libraries (R/sl_library.R) - SuperLearner arms only; plus T-learner (per-arm) DR outcome models |
| missing/binary | crossfit_rerun | 12600 | 12600 |    0 | 100.0 | complete | crossfitting strategy change; also the DGM was wrong three ways; plus scenario 2 (added later; was numbered 6 before the 2026-09-26 renumbering), rows 9901-12600; plus bug N (MNAR-Y truth), repaired at metrics time - no re-run; plus bug P and the risk-difference DGM - all 12,600 rows re-run, which also makes bug N moot; plus bug Q and the per-nuisance SuperLearner libraries (R/sl_library.R) - SuperLearner arms only; plus T-learner (per-arm) DR outcome models; plus the 2026-09-28 missing-data redesign - correlated covariates (X01-X03 auxiliary) and mechanisms MAR / MNAR-Y0 / MNAR-tau (missing/ADEMP.md) |
| missing/continuous | crossfit_rerun | 12600 | 12600 |    0 | 100.0 | complete | crossfitting strategy change; also bug F (dr_superlearner); plus bug O (continuous DGM) - all 12,600 rows re-run; plus scenario 2 (added later; was numbered 6 before the 2026-09-26 renumbering), rows 9901-12600; plus bug Q and the per-nuisance SuperLearner libraries (R/sl_library.R) - SuperLearner arms only; plus T-learner (per-arm) DR outcome models; plus the 2026-09-28 missing-data redesign - correlated covariates (X01-X03 auxiliary) and mechanisms MAR / MNAR-Y0 / MNAR-tau (missing/ADEMP.md) |
| crossfitting | dgm_rerun |   400 |   400 |    0 | 100.0 | complete | own comparison arms unchanged by the crossfitting change, but bug O changed the continuous DGM it runs on, and the 2026-09-26 renumbering its scenario ids (1/4/6/9 -> 1/4/6/8); old results archived; plus bug Q and the per-nuisance SuperLearner libraries (R/sl_library.R) - SuperLearner arms only |
| crossfitting/confidence_intervals | dgm_rerun |   150 |   150 |    0 | 100.0 | complete | pilot study, not part of the production rerun; re-run because bug O changed the continuous DGM and the 2026-09-26 renumbering its scenario ids (1/6/9 -> 1/4/8); old results archived |
| competing_risk/single_event | first_run |  3000 |  3000 |    0 | 100.0 | complete | new 2026-10-04: competing_risk/'s event 1 without event 2 - pseudo-value CF / DR / T, SuperLearner DR / T, RSF DR / T and the causal survival forest on the RMST, scenarios 1-3 (null, constant, heterogeneous), censoring on/off, 500 runs (competing_risk/single_event/README.md) |
| correlated/binary | first_run | 20800 | 20800 |    0 | 100.0 | complete | new 2026-10-01: binary/'s scenarios 1-4 with correlated covariates (copula, X01-X03 correlated too), at rho = 0 (the paired independent arm) and rho = 0.5, 500 runs (sample_size/correlated/README.md) |
| correlated/confidence_intervals/binary | first_run | 16000 | 16000 |    0 | 100.0 | complete | new 2026-10-04: confidence_intervals/binary/'s full CI_sf sweep on correlated/'s scenarios 1-4, rho = 0 and 0.5, n 500/1000, 100 runs (sample_size/correlated/confidence_intervals/README.md) |
| correlated/confidence_intervals/continuous | first_run | 16000 | 16000 |    0 | 100.0 | complete | new 2026-10-04: confidence_intervals/continuous/'s full CI_sf sweep on correlated/'s scenarios 1-4, rho = 0 and 0.5, n 500/1000, 100 runs (sample_size/correlated/confidence_intervals/README.md) |
| correlated/continuous | first_run | 20800 | 20800 |    0 | 100.0 | complete | new 2026-10-01: continuous/'s scenarios 1-4 with correlated covariates (copula, X01-X03 correlated too), at rho = 0 (the paired independent arm) and rho = 0.5, 500 runs (sample_size/correlated/README.md) |

## Legend

- **expected_jobs**: n_sims x number of parameter combinations in the study's grid
- **found_jobs**: res_sim_*.RDS files actually present under the study's results directory
- **status**: not_started (0 found), in_progress (0 < found < expected), complete (all found), blocked (currently fails to run, not scanned)

Generated by `Rscript check_all.R`. For resubmitting a specific study's
missing runs, use its own `<prefix>_check.R` -> `jobscripts/failed_ids.txt`
-> `qsub jobscripts/<prefix>_rerun.sh` loop; this script only reports.
