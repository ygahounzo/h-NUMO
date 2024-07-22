#!/bin/bash

# Compile the Fortran program (if not already compiled)
# if [ ! -f check ]; then
    gfortran -o check check.F90
# fi

# Run the Fortran executable
./check

cat check.log