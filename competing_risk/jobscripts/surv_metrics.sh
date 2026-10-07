#!/bin/bash
# Same mem rationale as surv_collect.sh - surv_metrics.R reads the same
# data-laden surv_all.RDS (~729MB serialized locally). ncpus = 2 for
# surv_metrics.R's `workers <- 2`; it strips each run to the truth and the
# framework estimates before shipping combos to the workers, so they hold far
# less than the parent.
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=32gb
#PBS -N surv_metrics
#PBS -j oe


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to the script directory
cd "${PBS_O_WORKDIR}/.."

Rscript surv_metrics.R
