#!/bin/bash
# resources updated after profiling runs, %max chunk updated since only 1 cpu 
# per run required (limit is 400 and this keeps space to start Rstudio OD)
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=1:ompthreads=1:mem=2gb
#PBS -J 1-2000%380
#PBS -N cf_1
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
Rscript cf_analysis.R "$PBS_ARRAY_INDEX"
