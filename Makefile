pr:
	mpicc knapsack_pruning_mpi.c -o knapsack_pruning_mpi
	./generator 100 20 | mpirun -n 6 ./knapsack_pruning_mpi
debug:
	mpicc -DDEBUG knap_MPI.c -o Knap_MPI
dp:
	mpicc knapsack_dp_mpi.c -o knapsack_dp_mpi
	./generator 100 20 | mpirun -n 6 ./knapsack_dp_mpi
compare:
	mpicc knapsack_dp_mpi.c -o knapsack_dp_mpi
	mpicc knapsack_pruning_mpi.c -o knapsack_pruning_mpi
	./generator 100 20 | mpirun -n 6 ./knapsack_dp_mpi
	./generator 100 20 | mpirun -n 6 ./knapsack_pruning_mpi
