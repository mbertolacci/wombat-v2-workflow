#!/bin/bash
#PBS -l wd

set -e
set -x

module load intel-compiler
module load netcdf/4.7.4
module load cdo

export GC_BIN=$NETCDF_ROOT/bin
export GC_INCLUDE=$NETCDF_ROOT/include/Intel
export GC_LIB=$NETCDF_ROOT/lib

BPCH_DIAG=n make -j24 mpbuild