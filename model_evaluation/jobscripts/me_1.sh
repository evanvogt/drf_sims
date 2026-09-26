#!/bin/bash
# PLACEHOLDER resources below - set by hand, never measured. Confirm it runs
# at all first:
#
#   Rscript model_evaluation/me_testing.R full
#
# The %N array throttle also still needs setting from the HPC queue's real
# memory/fair-share limits before the first real submission - each
# concurrent task starts its own H2O JVM cluster (mem="10G" heap), which
# rules out anything resembling continuous's %190 or crossfitting's %380.
# See README.md.
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=10gb
#PBS -J 1-360%4
#PBS -N me_1
#PBS -o logs_1/
#PBS -e logs_1/


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters. Trailing args are workers/n_cores, set by
# hand to match ncpus above - change them together.
Rscript me_analysis.R "$PBS_ARRAY_INDEX" 2 2
