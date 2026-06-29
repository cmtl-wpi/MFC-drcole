# nighthawk GPU build/run env (NVHPC 24.11 + bundled hpcx MPI, single V100).
# Source before every ./mfc.sh build/run with --gpu acc. See memory/gpu-build-on-nighthawk.md
export NVHPC=/opt/nvidia/hpc_sdk/Linux_x86_64/24.11
export PATH=$NVHPC/comm_libs/mpi/bin:$NVHPC/compilers/bin:$PATH
export LD_LIBRARY_PATH=$NVHPC/comm_libs/12.6/hpcx/hpcx-2.20/ompi/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/lib64:$NVHPC/cuda/12.6/targets/x86_64-linux/lib:$NVHPC/cuda/12.6/lib64:$NVHPC/math_libs/12.6/targets/x86_64-linux/lib:$LD_LIBRARY_PATH
export CC=nvc CXX=nvc++ FC=nvfortran OMPI_CC=nvc OMPI_CXX=nvc++ OMPI_FC=nvfortran
# Runtime: don't bind to a single core; oversubscribe the single GPU for n>1
export OMPI_MCA_hwloc_base_binding_policy=none
export OMPI_MCA_rmaps_base_oversubscribe=1
