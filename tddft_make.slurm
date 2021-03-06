#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=75:00:00
#SBATCH --error=%J.err
#SBATCH --output=%J.out
#SBATCH --job-name=TDDFT

ml ic

make

./a.out 0.5 10.0 0.0
#./a.out 0.5 10.0 1.0
#./a.out 0.5 10.0 2.0
./a.out 0.5 10.0 3.0


./a.out 1.0 10.0 0.0
#./a.out 1.0 10.0 1.0
#./a.out 1.0 10.0 2.0
./a.out 1.0 10.0 3.0

./a.out 1.5 10.0 0.0
#./a.out 1.5 10.0 1.0
#./a.out 1.5 10.0 2.0
./a.out 1.5 10.0 3.0


make clean

#rm  3*




