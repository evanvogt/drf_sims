#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=4:ompthreads=4:mem=8gb
#PBS -N cf_collect


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to the script directory
cd "${PBS_O_WORKDIR}/.."

# reads every per-run file and computes the metrics in one streaming pass.
# modest memory because the per-run files carry no data or nuisance matrices -
# roughly 0.4 mb each, against the 15gb the other studies' collect jobs need.
Rscript cf_collect.R
