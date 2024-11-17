#!/bin/bash
###
###
#SBATCH -J bump             # job name
#SBATCH -o log.o%j         # output and error file name (%j expands to jobID)
#SBATCH -n 40                # total number of tasks requested
#SBATCH -N 1                # number of nodes you want to run on
#SBATCH --cpus-per-task 1
#SBATCH -p coaps_q          # queue (partition)
#SBATCH -t 01:00:00         # run time (hh:mm:ss)
#SBATCH --exclusive

cp ~/h-NUMO/bin/numo3d .

mpirun -np 1 ./numo3d
