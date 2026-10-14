#!/bin/bash

#SBATCH --job-name="test"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -t 03:00
##SBATCH --gres=gpu:1
#SBATCH --gres=gpu:ada4500:1
#SBATCH -A backfill2
#SBATCH -o out.log
#SBATCH -e err.log
#SBATCH --mail-type=NONE

# module load nvhpc/26.5
# NVHPC_MPI=/gpfs/research/software/NVHPC/26.5/Linux_x86_64/26.5/comm_libs/mpi/bin
# export PATH="$NVHPC_MPI:$PATH"

NV=/gpfs/research/software/NVHPC/26.5/Linux_x86_64/26.5
export PATH=$NV/comm_libs/mpi/bin:$NV/compilers/bin:$PATH
export LD_LIBRARY_PATH=$NV/compilers/lib:$NV/cuda/lib64:$LD_LIBRARY_PATH

cp ../bin/numo3d .

# srun -n 2 ./numo3d
mpirun -np 1 ./numo3d

