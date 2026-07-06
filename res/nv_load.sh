NVARCH=Linux_x86_64
NVCOMPILERS=/opt/nvidia/hpc_sdk

export NVARCH
export NVCOMPILERS

export PATH=$NVCOMPILERS/$NVARCH/26.3/compilers/bin:$PATH
export PATH=$NVCOMPILERS/$NVARCH/26.3/comm_libs/mpi/bin:$PATH
export MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/26.3/compilers/man


# CaNS cuDecomp runtime library
CANS_ROOT=$HOME/codes/CaNS-IBM
export CANS_ROOT
export LD_LIBRARY_PATH=$CANS_ROOT/dependencies/cuDecomp/build/lib:$LD_LIBRARY_PATH
