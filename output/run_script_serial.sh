#!/bin/bash

cp ../bin/numo3d .

mpirun -np 4 ./numo3d
# ./numo3d
