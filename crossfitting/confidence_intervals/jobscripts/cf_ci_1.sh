#!/bin/bash
# Pilot: half-sample bootstrap CIs for the 11 RF/CF arms.
# 3 scenarios x 50 runs = 150 array jobs. Resources below are hand-set
# placeholders - see crossfitting/README.md's "Sizing the CI pilot's array job"
# before submitting the full array. The bootstrap refits (~200 x V forests per
# crossfit arm x 4 crossfit arms) dominate cost here, well beyond what cf_1.sh's
# point-estimate-only timings would suggest.
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=5gb
#PBS -J 1-150
#PBS -N cf_ci_1
#PBS -o logs_ci_1/
#PBS -e logs_ci_1/


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters. Trailing args are CI_boot/workers/grf_threads,
# set by hand to match ncpus above - change them together.
Rscript cf_ci_analysis.R "$PBS_ARRAY_INDEX" 200 2 1
