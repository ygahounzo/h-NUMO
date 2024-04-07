#!/bin/bash
###
###
#SBATCH --time=24:00:00
#SBATCH -n 48
#SBATCH --partition=bsudfq
#SBATCH --job-name=N4Ne25free
#SBATCH --output=output.o%j
#SBATCH -N 1
#SBATCH --exclusive

source /etc/profile
source ~/.bashrc
source ~/.bash_profile

cp ~/h-NUMO/bin/numo3d .

mpirun -np 48 ./numo3d