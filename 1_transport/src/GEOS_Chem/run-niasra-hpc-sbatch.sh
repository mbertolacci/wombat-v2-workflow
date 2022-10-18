#!/bin/bash
#SBATCH --hint=nomultithread
#SBATCH --output=output/run.log

set -aex

export LD_LIBRARY_PATH="/opt/intel/oneapi/vpl/2022.0.0/lib:/opt/intel/oneapi/tbb/2021.5.0/env/../lib/intel64/gcc4.8:/opt/intel/oneapi/mpi/2021.5.0//libfabric/lib:/opt/intel/oneapi/mpi/2021.5.0//lib/release:/opt/intel/oneapi/mpi/2021.5.0//lib:/opt/intel/oneapi/mkl/2022.0.1/lib/intel64:/opt/intel/oneapi/itac/2021.5.0/slib:/opt/intel/oneapi/ipp/2021.5.1/lib/intel64:/opt/intel/oneapi/ippcp/2021.5.0/lib/intel64:/opt/intel/oneapi/ipp/2021.5.1/lib/intel64:/opt/intel/oneapi/dnnl/2022.0.1/cpu_dpcpp_gpu_dpcpp/lib:/opt/intel/oneapi/debugger/2021.5.0/gdb/intel64/lib:/opt/intel/oneapi/debugger/2021.5.0/libipt/intel64/lib:/opt/intel/oneapi/debugger/2021.5.0/dep/lib:/opt/intel/oneapi/dal/2021.5.1/lib/intel64:/opt/intel/oneapi/compiler/2022.0.1/linux/lib:/opt/intel/oneapi/compiler/2022.0.1/linux/lib/x64:/opt/intel/oneapi/compiler/2022.0.1/linux/lib/oclfpga/host/linux64/lib:/opt/intel/oneapi/compiler/2022.0.1/linux/compiler/lib/intel64_lin:/opt/intel/oneapi/ccl/2021.5.0/lib/cpu_gpu_dpcpp"

ulimit -s unlimited
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_STACKSIZE=500m
srun -c $OMP_NUM_THREADS ./geos.mp > output/geos.mp.log

source /etc/profile.d/conda.sh
conda activate ../../../../.conda_env
srun -c $OMP_NUM_THREADS bash postprocess-run.sh
