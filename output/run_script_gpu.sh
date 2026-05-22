#!/bin/bash

#SBATCH --job-name="test"
#SBATCH --ntasks=1
#SBATCH -t 1:00
#SBATCH --gres=gpu:1
#SBATCH -A backfill2
#SBATCH -o out.log
#SBATCH -e err.log
#SBATCH --mail-type=NONE



module load nvhpc/23.11
module load openmpi

cp ../bin/numo3d .

# echo "PATH=$PATH"
# echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

# which nvfortran
# which mpirun
# which orted || true
# which pmix_info || true

# mpirun --version
# ompi_info | head -20

srun -n 2 ./numo3d
