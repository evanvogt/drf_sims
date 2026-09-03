#!/bin/bash
# Same mem rationale as surv_collect.sh - surv_metrics.R reads the same
# data-laden surv_all.RDS (~729MB serialized locally) and unnests it further.
# Single core: unlike surv_collect.R's get_results(workers = 2), this script
# does no parallel work.
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:ompthreads=1:mem=32gb
#PBS -N surv_metrics


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to the script directory
cd "${PBS_O_WORKDIR}/.."

Rscript surv_metrics.R
