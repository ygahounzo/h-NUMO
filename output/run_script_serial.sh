#!/bin/bash

cp ../bin/numo3d .

mpirun -np 1 ./numo3d
# ./numo3d
