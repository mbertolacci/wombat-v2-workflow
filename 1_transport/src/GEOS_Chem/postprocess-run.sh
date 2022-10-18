#!/bin/bash

set -ex

mv HEMCO_restart.*.nc output

for month in %1$s; do
    rm output/HEMCO_diagnostics.${month}*
done

if ls output | grep -q HEMCO_diagnostics; then
    cdo -L -z zip_6 -W \
        -splityearmon \
        -cat output/HEMCO_diagnostics\*nc \
        output/HEMCO_diagnostics.
fi

find output/ \
    -name "HEMCO_diagnostics.????????????.nc" \
    -delete

for field in %2$s; do
    for month in %3$s; do
        find output/ \
            -name "GEOSChem.${field}Hourly.${month}??_0000z.nc4" \
            -delete
    done
done

gzip -f output/HEMCO.log output/geos.mp.log

ln -sf output/GEOSChem.Restart.* .
