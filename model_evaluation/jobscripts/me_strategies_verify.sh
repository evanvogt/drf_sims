#!/bin/bash
# Proves me_strategies.R's pass-through was inert: the 9 candidate fits, the
# data, the truth, fold_info and the `whole`/`cv_indep` arms must be
# bit-identical to the source tree in every run. Pure I/O over ~716 files, no
# model fitting - minutes, not hours.
#
# Run after me_strategies.sh finishes and before me_collect.sh.
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=16gb
#PBS -N me_strat_verify


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to the script directory
cd "${PBS_O_WORKDIR}/.."

Rscript me_strategies_verify.R
