#!/bin/bash
#PBS -l walltime=01:00:00  
#PBS -l select=1:ncpus=2:ompthreads=2:mem=100gb
#PBS -N ci_miss_metrics

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script
Rscript cts_miss_metrics.R