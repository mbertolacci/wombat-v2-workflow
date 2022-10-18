#!/bin/bash

set -e
set -x
set -o pipefail

module load \
    R \
    intel-compiler \
    netcdf/4.7.4 \
    intel-mkl \
    udunits \
    gdal \
    proj \
    geos

export WOMBAT_LOG_LEVEL=trace
export WOMBAT_MAX_WORKERS=$PBS_NCPUS
export UTILS_PARTIAL={utils_partial}

Rscript 2_matching/src/match.R \
    --oco2-soundings {oco2_soundings} \
    --tccon-sounding-directory {tccon_sounding_directory} \
    --obspack-directory {obspack_directory} \
    --lauder {lauder} \
    --input-run {input_run} \
    --meteorology-run {meteorology_run} \
    --output {output}
