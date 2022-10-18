#!/bin/bash
#PBS -l wd

set -e
set -x
set -o pipefail

module load intel-compiler
module load netcdf/4.7.4
module load cdo

ulimit -s unlimited
export OMP_NUM_THREADS=$PBS_NCPUS
export OMP_STACKSIZE=500m

# Set up output directory on jobfs
rm -rf output
ln -s $PBS_JOBFS output
./geos.mp | tee output/geos.mp.log

bash postprocess-run.sh

# Move outputs out of jobfs
mkdir -p output2
rsync -a output/ output2/
rm output
mv output2 output
