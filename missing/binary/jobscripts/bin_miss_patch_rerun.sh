#!/bin/bash
# Re-run the combinations bin_miss_patch_check.R found unaccounted for.
#
# Indices come from failed_patch_ids.txt, one per line, exactly as
# bin_miss_rerun.sh reads failed_ids.txt - but note the two lists index
# DIFFERENT things. failed_ids.txt holds rows of study$grid (one simulation
# run); this holds rows of combos(study) (one parameter combination and its 100
# files), which is what bin_miss_patch.R takes.
#
# Set -J to 1-<number of lines in failed_patch_ids.txt> before submitting;
# bin_miss_patch_check.R prints the number. There is no update_rerun_script()
# for this job - that function rewrites <prefix>_rerun.sh against
# <prefix>_1.sh, and neither is this pair.
#
# Safe to submit over combinations that are already correct: patch_status_of()
# recognises an already-patched file and the element just rewrites its manifest.
#
# Three deliberate differences from bin_miss_patch.sh, which is how the original
# ten elements were lost:
#   - a %N throttle. bin_miss_patch.sh has none, so all 99 elements ran at once,
#     each doing 100 readRDS + 100 saveRDS against the shared filesystem.
#   - double the walltime, following the time_factor = 2 convention in
#     R/pipeline.R's update_rerun_script().
#   - -j oe, so stderr is merged into the .o file. The HPC is not returning .e
#     files at the moment, which is why the original failure left no evidence.
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=1:ompthreads=1:mem=8gb
#PBS -J 1-10%4
#PBS -N bin_miss_patch_rerun
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
Rscript bin_miss_patch.R "$comboid"
