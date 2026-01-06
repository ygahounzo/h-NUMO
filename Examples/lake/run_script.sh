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

# NB: Adjust number of processes (-np) according to the number of cores allocated by SLURM on your system.

cp $(pwd)/../../bin/numo3d .

# mpirun -np 36 ./numo3d
./numo3d
