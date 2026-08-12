#!/bin/bash
# the full 9900-row grid from bin_miss_config.R. This used to be 1-1100, the
# size AFTER bin_miss_analysis.R filtered to the complete_data arm - the filter
# that renumbered every index (bug D). To run one arm, use
#   grid_indices(study, method = "complete_data")
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=1:ompthreads=1:mem=2gb
#PBS -J 1-9900%190
#PBS -N bin_miss_1
#PBS -o logs_1/
#PBS -e logs_1/

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript bin_miss_analysis.R "$PBS_ARRAY_INDEX" 1 1
