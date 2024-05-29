#!/bin/bash
###
###
#SBATCH --time=01:00:00
#SBATCH -n 36
#SBATCH --partition=bsudfq
#SBATCH --job-name=bump
#SBATCH --output=output.o%j
#SBATCH -N 1
#SBATCH --exclusive

cp ~/h-NUMO/bin/numo3d .

mpirun -np 36 ./numo3d
