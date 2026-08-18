#!/bin/bash
# Re-run the combinations cts_miss_patch_check.R found unaccounted for.
# The twin of missing/binary/jobscripts/bin_miss_patch_rerun.sh, which carries
# the full explanation of why this pair exists.
#
# This study's patch completed (check_all.R has it at 8,800/8,800), so this is
# here for symmetry and for the next time - not because anything is outstanding.
#
# Indices come from failed_patch_ids.txt and are rows of combos(study) - one
# parameter combination and its 100 files - NOT rows of study$grid, which is
# what failed_ids.txt and cts_miss_rerun.sh hold. Set -J to
# 1-<number of lines in failed_patch_ids.txt> before submitting;
# cts_miss_patch_check.R prints the number.
#
# Safe to submit over combinations that are already correct: the patch is
# idempotent, so an element that only lost its manifest just rewrites it.
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=1:ompthreads=1:mem=8gb
#PBS -J 1-1%4
#PBS -N cts_miss_patch_rerun
#PBS -o logs_patch/
#PBS -e logs_patch/
#PBS -j oe

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# failed combination ids
cd "${PBS_O_WORKDIR}"
comboid=$(sed -n "${PBS_ARRAY_INDEX}p" failed_patch_ids.txt)

echo "re-patching combination: ${comboid}"

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript cts_miss_patch.R "$comboid"
