#!/bin/bash
#PBS -l walltime=3:00:00
#PBS -l select=1:ncpus=10:ompthreads=10:mem=20gb
#PBS -J 1-249%10
#PBS -N cmc_rerun
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

# Run Rscript with parameters. Trailing args are workers (imputations fitted
# in parallel)/grf_threads: one worker per cpu, one grf thread each, so
# workers should match ncpus above - change them together.
Rscript cts_miss_ci_analysis.R "$jobid" 10 1
