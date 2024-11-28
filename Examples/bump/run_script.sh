#!/bin/bash
###
###
#SBATCH --time=01:00:00
#SBATCH -n 40
#SBATCH --partition=coaps_q
#SBATCH --job-name=bump
#SBATCH --output=output_log.o%j
#SBATCH -N 1
#SBATCH --exclusive

cp $(pwd)/../../bin/numo3d .

mpirun -np 36 ./numo3d
