#!/bin/bash
#PBS -l walltime=02:00:00  
#PBS -l select=1:ncpus=2:ompthreads=2:mem=10gb
#PBS -N cts_rerun
#PBS -o logs_rerun/
#PBS -e logs_rerun/


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript cts_analysis.R "1402"