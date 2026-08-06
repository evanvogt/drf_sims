#!/bin/bash
# resources updated after profiling runs, %max chunk updated since only 1 cpu
# per run required (limit is 400 and this keeps space to start Rstudio OD)
# NOTE: the trailing workers/grf_threads args below (2 1) were never re-synced to
# the last profiling run that set ncpus=1 above - the mismatch cf_profile_summary.R
# now avoids by writing both together. Re-run the profiling sweep before the next
# production submission to fix this pair.
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

# Run R script with parameters. Trailing args are workers/grf_threads - placeholders
# here, overwritten by cf_profile_summary.R alongside the #PBS -l lines above.
Rscript cf_analysis.R "$PBS_ARRAY_INDEX" 2 1
