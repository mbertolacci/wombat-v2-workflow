#!/bin/bash

set -aex

input_file=$1
output_file=$2

tmp_file=$(mktemp)
tmp_file2=$(mktemp)

cdo -w -z zip_6 \
    -shifttime,-14days \
    -setattribute,fgco2_raw@units=kg/m2/s \
    -divdpy \
    -divc,86400 \
    -mulc,0.04401 \
    -setmisstoc,0 \
    -select,name=fgco2_raw \
    $input_file \
    $tmp_file

ncwa --overwrite \
    -a bnds \
    $tmp_file \
    $tmp_file2

ncks -C -O -x \
    -v time_bnds,lon_bnds,lat_bnds \
    $tmp_file2 \
    $output_file

ncatted \
    -a Title,global,c,c,"Ocean CO2 flux created from NCEI accession doi.org/10.7289/V5Z899N6" \
    -a bounds,time,d,, \
    -a bounds,lat,d,, \
    -a bounds,lon,d,, \
    $output_file

ncatted \
    -a _FillValue,fgco2_raw,d,, \
    -a missing_value,fgco2_raw,d,, \
    -a FillValue,fgco2_raw,d,, \
    $output_file
