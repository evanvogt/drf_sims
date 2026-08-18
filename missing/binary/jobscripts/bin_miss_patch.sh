#!/bin/bash
# Back-fill the dr_random_forest HTE tests into the finished results.
# One-off repair, not part of the simulation pipeline - see R/patch_hte_tests.R.
#
# The array index is a row of combos(study), i.e. one parameter combination and
# its 100 res_sim_*.RDS files - NOT a row of study$grid, which is what
# bin_miss_1.sh indexes. There are 99 combinations:
#   nrow(combos(study))  after sourcing bin_miss_config.R
#
# Each element reads 100 files, runs one OLS and two asymptotic permutation
# tests per file at n = 500, and rewrites them; ~30-60 s in total. The walltime
# below is deliberately far above that - it costs nothing in the queue and this
# job rewrites irreplaceable results, so it must never be killed part-way.
# (It could not corrupt one if it were: the write is saveRDS-to-.tmp then
# rename. But a half-done combination is still a manifest to reconcile.)
#
# Safe to re-run: already-patched files are detected and skipped.
# Dry run first - `Rscript bin_miss_patch.R dry` on the login node is cheap.
#PBS -l walltime=00:20:00
#PBS -l select=1:ncpus=1:ompthreads=1:mem=4gb
#PBS -J 1-99
#PBS -N bin_miss_patch
#PBS -o logs_patch/
#PBS -e logs_patch/

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript bin_miss_patch.R "$PBS_ARRAY_INDEX"
