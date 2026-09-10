#!/bin/bash
#PBS -l walltime=00:15:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=3gb
#PBS -J 2000-7000%190
#PBS -N surv_extra
#PBS -o logs_extra/
#PBS -e logs_extra/


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript surv_analysis.R "$PBS_ARRAY_INDEX"