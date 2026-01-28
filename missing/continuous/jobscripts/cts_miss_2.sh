#!/bin/bash
#PBS -l walltime=02:00:00  
#PBS -l select=1:ncpus=2:ompthreads=2:mem=10gb
#PBS -J 10001-12200%190
#PBS -N cts_miss_2
#PBS -o logs_2/
#PBS -e logs_2/
#PBS -W depend=afterany:1618334[]

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env


# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript cts_miss_analysis.R "$PBS_ARRAY_INDEX"
