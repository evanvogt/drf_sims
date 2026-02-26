#!/bin/bash
#PBS -l walltime=02:00:00  
#PBS -l select=1:ncpus=2:ompthreads=2:mem=200gb
#PBS -N ci_cts_metrics
#PBS -W depend=afterok:1612866

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script
Rscript cts_ci_metrics.R