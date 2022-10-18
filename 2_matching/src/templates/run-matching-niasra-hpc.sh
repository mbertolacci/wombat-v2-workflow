#!/bin/bash
#SBATCH --hint=nomultithread
#SBATCH --output={output}/run-matching.log

set -e
set -x
set -o pipefail

source /etc/profile.d/conda.sh
conda activate ./.conda_env

export WOMBAT_LOG_LEVEL=trace
export WOMBAT_MAX_WORKERS=$SLURM_CPUS_PER_TASK
export UTILS_PARTIAL={utils_partial}

Rscript 2_matching/src/match.R \
    --oco2-soundings {oco2_soundings} \
    --tccon-sounding-directory {tccon_sounding_directory} \
    --obspack-directory {obspack_directory} \
    --lauder {lauder} \
    --input-run {input_run} \
    --meteorology-run {meteorology_run} \
    --output {output}
