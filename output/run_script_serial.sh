#!/bin/bash

cp ../bin/numo3d .

mpirun -np 8 ./numo3d
# ./numo3d
