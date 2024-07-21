#!/bin/bash
###
###
#SBATCH --time=01:00:00
#SBATCH -n 48
#SBATCH --partition=bsudfq
#SBATCH --job-name=sc
#SBATCH --output=output.o%j
#SBATCH -N 1
#SBATCH --exclusive

cp ~/NUMO_MLSWE/MultiLayers/h-NUMO/bin/numo3d .

mpirun -np 1 ./numo3d
