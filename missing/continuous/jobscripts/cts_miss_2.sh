#!/bin/bash
#PBS -l walltime=02:00:00  
#PBS -l select=1:ncpus=2:ompthreads=2:mem=10gb
# 5001-9900 completes the 9900-row grid in cts_miss_config.R.
# This used to be 10001-12200: indices that never existed in any version of the
# grid. See the note in cts_miss_1.sh.
#PBS -J 5001-9900%190
#PBS -N cts_miss_2
#PBS -o logs_2/
#PBS -e logs_2/
# NOTE: this used to carry `#PBS -W depend=afterany:1618334[]`, a dependency on
# a long-finished job. Left in, it holds the array forever. Re-add a dependency
# only with the CURRENT job id if you want the two halves to run in sequence.

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env


# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript cts_miss_analysis.R "$PBS_ARRAY_INDEX"
