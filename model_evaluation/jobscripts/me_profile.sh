#!/bin/bash
# Deliberately over-provisioned so nothing fails for the wrong reason.
# ncpus=8 covers the largest configuration in the sweep (n_cores=5, workers=4
# - these parallelise sequential phases, not concurrent ones, so the ceiling
# is max(workers, n_cores) + headroom, not their sum). mem=24gb covers H2O's
# own "10G" heap plus R/XGBoost/SuperLearner working memory above it.
#
# Unlike crossfitting/jobscripts/cf_profile.sh's unthrottled -J 1-36 (cheap
# RF/causal-forest fits, no external process), each cell here starts its own
# H2O JVM cluster - throttled to %4 so the profiling sweep itself doesn't hit
# the same per-task-H2O-cluster memory pressure the real array job needs to
# avoid. See me_profile.R for what each cell measures.
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=8:ompthreads=8:mem=24gb
#PBS -J 1-16%4
#PBS -N me_profile
#PBS -o logs_profile/
#PBS -e logs_profile/


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

Rscript me_profile.R "$PBS_ARRAY_INDEX"
