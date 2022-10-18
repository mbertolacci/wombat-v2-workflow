#!/bin/bash
#SBATCH --hint=nomultithread
#SBATCH --output={output}/run-aggregate-flux.log

set -e
set -x
set -o pipefail

source /etc/profile.d/conda.sh
conda activate ./.conda_env

cdo -v -W -L -z zip_6 -monavg -vertsum \
    -cat {input_run}\*/output/HEMCO_diagnostics.\*.nc \
    {output}/monthly-fluxes.nc
