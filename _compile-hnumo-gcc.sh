#!/bin/bash
. /etc/bashrc
module load compile/gcc/7.2.0
module load openmpi/2.0.4
module load netcdf-fortran/4.4.4
module load lib/atlas+lapack/3.10.2+3.5.0
make hnumo-hamming-p4est-gcc7
