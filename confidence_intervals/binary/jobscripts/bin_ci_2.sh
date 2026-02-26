#!/bin/bash
#PBS -l walltime=02:00:00  
#PBS -l select=1:ncpus=2:ompthreads=2:mem=5gb
#PBS -J 10001-20000%190
#PBS -N ci_bin_2
#PBS -o logs_2/
#PBS -e logs_2/
#PBS -W depend=afterany:1626671[]

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript bin_ci_analysis.R "$PBS_ARRAY_INDEX"