#To run, just type "compile-hnumo-gcc.sh"
#!/bin/bash
#DIR is where the script is located (at $NUMO_DIR but command below allows it to be anywhere)                                       
DIR=$(dirname "$0")
srun -n 1 --time=00:59:00 bash --noprofile --norc $DIR/_compile-hnumo-gcc.sh