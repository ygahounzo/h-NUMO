#!/bin/bash

#SBATCH --job-name="test"
#SBATCH --ntasks=1
#SBATCH -t 03:00
#SBATCH --gres=gpu:1
#SBATCH -A backfill2
#SBATCH -o out.log
#SBATCH -e err.log
#SBATCH --mail-type=NONE

module load nvhpc/23.11
module load openmpi

cp ../bin/numo3d .

srun -n 1 ./numo3d
