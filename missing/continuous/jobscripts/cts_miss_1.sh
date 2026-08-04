#!/bin/bash
#PBS -l walltime=03:00:00  
#PBS -l select=1:ncpus=2:ompthreads=2:mem=10gb
# 1-5000 of the 9900-row grid in cts_miss_config.R; cts_miss_2.sh does the rest.
# This used to be 1-1100, which was the size of the grid AFTER
# cts_miss_analysis.R filtered it to the complete_data arm - the filter that
# renumbered every index (bug D). The grid is no longer filtered, so the range
# is the whole design. To run one arm only, use
#   grid_indices(study, method = "complete_data")
# and submit those indices.
#PBS -J 1-5000%190
#PBS -N cts_miss_1
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
Rscript cts_miss_analysis.R "${PBS_ARRAY_INDEX}"
