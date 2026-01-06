#!/bin/bash
###
###
#SBATCH --time=12:00:00
#SBATCH -n 40
#SBATCH --partition=coaps_q
#SBATCH --job-name=gyre
#SBATCH --output=output_log.o%j
#SBATCH -N 1
#SBATCH --exclusive

cp $(pwd)/../../bin/numo3d .

mpirun -np 40 ./numo3d
