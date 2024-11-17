#!/bin/bash

# Remove the compiled Fortran program if already existed
if [ -f check ]; then
    rm check
fi

./numo3d
# Compile the Fortran program
gfortran -o check check.F90

# Run the Fortran executable
./check

cat output.log