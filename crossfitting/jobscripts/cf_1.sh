#!/bin/bash
#PBS -l walltime=00:40:00
#PBS -l select=1:ncpus=1:ompthreads=1:mem=2gb
#PBS -J 1-400
#PBS -N cf_1
#PBS -o logs_1/
#PBS -j oe


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters. Trailing args are workers/grf_threads, set by
# hand to match ncpus above - change them together.
Rscript cf_analysis.R "$PBS_ARRAY_INDEX" 1 1
