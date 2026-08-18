#!/bin/bash
# Back-fill the dr_random_forest HTE tests into the finished results.
# One-off repair, not part of the simulation pipeline - see R/patch_hte_tests.R.
# As missing/binary/jobscripts/bin_miss_patch.sh; that script's header has the
# full reasoning.
#
# The array index is a row of combos(study) - one parameter combination and its
# 100 res_sim_*.RDS files - NOT a row of study$grid, which is what cts_miss_1.sh
# indexes. There are 99 combinations.
#
# Submit this only once cts_miss_rerun.sh has finished and cts_miss_check.R
# reports 9,900/9,900, so the patch is one clean pass over a complete study.
#
# Safe to re-run: already-patched files are detected and skipped.
# Dry run first - `Rscript cts_miss_patch.R dry` on the login node is cheap.
#PBS -l walltime=00:20:00
#PBS -l select=1:ncpus=1:ompthreads=1:mem=4gb
#PBS -J 1-99%20
#PBS -N cts_miss_patch
#PBS -o logs_patch/
#PBS -e logs_patch/
#PBS -j oe

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript cts_miss_patch.R "$PBS_ARRAY_INDEX"
