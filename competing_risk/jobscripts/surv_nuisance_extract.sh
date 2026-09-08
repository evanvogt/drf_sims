#!/bin/bash
# ncpus=7/ompthreads=7 must match surv_nuisance_extract.R's `workers <- 7`
# (future::multisession, one worker per (scenario, censoring) combo - 14 of
# them, so this could go up to 14 given more cores). mem=48gb: same baseline
# as surv_collect.sh's 32gb for holding surv_all.RDS unserialized (~729MB
# serialized locally, larger in memory), plus headroom for 7 worker
# processes each holding their own ~1/14th slice plus working tibbles - a
# reasoned estimate, not a profiled one (this script is new and untested on
# the cluster). walltime=02:00:00 similarly: surv_metrics.R's single-core
# equivalent pass over the same 7,000 runs gets 1h, but this script does more
# per-run work (building and pivoting a long tibble per nuisance arm x
# estimand cell) before the 7x parallel speedup; check actual job stats after
# the first submission and tighten both.
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
