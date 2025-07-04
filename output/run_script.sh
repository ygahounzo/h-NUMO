#!/bin/bash
###
###
#SBATCH -J bump             # job name
#SBATCH -o out.log         # output and error file name (%j expands to jobID)
#SBATCH -n 40                # total number of tasks requested
##SBATCH -N 1                # number of nodes you want to run on
##SBATCH --cpus-per-task 1
#SBATCH -p coaps_q          # queue (partition)
#SBATCH -t 00:10:00         # run time (hh:mm:ss)
#SBATCH --exclusive

cp ../bin/numo3d .

srun -n 25 ./numo3d
