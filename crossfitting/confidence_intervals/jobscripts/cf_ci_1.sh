#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=5gb
#PBS -J 1-150
#PBS -N cf_ci_1
#PBS -o logs_ci_1/
#PBS -j oe


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters. Trailing args are CI_boot/workers/grf_threads,
# set by hand to match ncpus above - change them together.
Rscript cf_ci_analysis.R "$PBS_ARRAY_INDEX" 200 2 1
