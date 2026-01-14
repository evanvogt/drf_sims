#!/bin/bash
#PBS -l walltime=01:00:00  
#PBS -l select=1:ncpus=2:ompthreads=2:mem=5gb
#PBS -J 1-9
#PBS -N cts_miss_rerun
#PBS -o logs_rerun/
#PBS -e logs_rerun/

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# failed job ids
cd "${PBS_O_WORKDIR}"
jobid=$(sed -n "${PBS_ARRAY_INDEX}p" failed_ids.txt)

echo "rerunning index: $jobid"

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript cts_miss_analysis.R "$jobid"
