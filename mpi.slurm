#!/bin/bash
#SBATCH --time=0-01:30:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
module load OpenMPI/2.0.0-GCC-6.2.0
mpicc knapsack_pruning_mpi.c -lm -o main
mpiexec -n 4 ./main < test5000300.txt