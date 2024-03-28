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

source /etc/profile
source ~/.bashrc
source ~/.bash_profile

cp ~/numo/bin/numa3d .

mpirun -np 48 ./numa3d
