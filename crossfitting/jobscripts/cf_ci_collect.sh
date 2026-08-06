#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=4:ompthreads=4:mem=4gb
#PBS -N cf_ci_collect


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to the script directory
cd "${PBS_O_WORKDIR}/.."

# reads every per-run file and computes the metrics in one streaming pass.
# only 150 replicates x 5 arms here, so 4gb is generous headroom.
Rscript cf_ci_collect.R
