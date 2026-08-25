#!/bin/bash
# The FULL verification suite, on the cluster. me_testing.R's `full` mode is
# the slow half (real XGBoost-CV and H2O AutoML fits) and SL2 fails on the
# local Windows machine for a documented package-version reason, so `full` is
# not a useful local signal - run it here instead, and run only the fast form
# locally. See README.md.
#
# Submit this BEFORE me_1.sh / me_strategies.sh. As well as checking the code,
# it settles two things about the environment that nothing else does: whether
# sim-env really carries h2o/xgboost/caret with a working Java runtime, and
# whether SL2's local failure reproduces here (it should not).
#
# Not an array job, and sized for a single H2O JVM (mem=10G heap) plus the
# candidate fits at n=250.
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=4:ompthreads=2:mem=12gb
#PBS -N me_testing
#PBS -o logs_testing/
#PBS -e logs_testing/


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to the script directory
cd "${PBS_O_WORKDIR}/.."

Rscript me_testing.R full
