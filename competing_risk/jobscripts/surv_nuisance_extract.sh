#!/bin/bash
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=7:ompthreads=7:mem=48gb
#PBS -N surv_nuisance_extract


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to the script directory
cd "${PBS_O_WORKDIR}/.."

Rscript surv_nuisance_extract.R
