#!/bin/bash
# One grid row per subjob, from failed_ids_oracle.txt. study.R check (on
# snr_oracle_config.R) rewrites -J to match it, and lifts the resources above
# snr_oracle_1.sh's - far more than one row needs, which is harmless.
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=4gb
#PBS -J 1-2%100
#PBS -N snr_oracle_rerun
#PBS -o logs_oracle_rerun/
#PBS -j oe

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

jobid=$(sed -n "${PBS_ARRAY_INDEX}p" "${PBS_O_WORKDIR}/failed_ids_oracle.txt")

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters: one row, grf threads
Rscript snr_oracle_analysis.R "$jobid" 1 2
