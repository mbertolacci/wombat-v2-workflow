#!/bin/bash

set -e
set -x
set -o pipefail

module load cdo

cdo -v -W -L -z zip_6 -monavg -vertsum \
    -cat {input_run}\*/output/HEMCO_diagnostics.\*.nc \
    {output}/monthly-fluxes.nc
