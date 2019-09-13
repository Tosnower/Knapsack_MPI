pr:
	mpicc knapsack_pruning_mpi.c -o knapsack_pruning_mpi
	./generator 100 20 | mpirun -n 6 ./knapsack_pruning_mpi
prc:
	mpicc knapsack_pruning_mpi.c -o knapsack_pruning_mpi
	cat test30.txt | mpirun -n 6 ./knapsack_pruning_mpi

debug:
	mpicc -DDEBUG knap_MPI.c -o Knap_MPI
dp:
	mpicc knapsack_dp_mpi.c -o knapsack_dp_mpi
	./generator 100 20 | mpirun -n 6 ./knapsack_dp_mpi
com:
	mpicc knapsack_dp_mpi.c -o knapsack_dp_mpi
	mpicc knapsack_pruning_mpi.c -o knapsack_pruning_mpi
	./generator 1000 200 | mpirun -n 6 ./knapsack_dp_mpi
	./generator 1000 200 | mpirun -n 6 ./knapsack_pruning_mpi
comf:
		mpicc knapsack_dp_mpi.c -o knapsack_dp_mpi
		mpicc knapsack_pruning_mpi.c -o knapsack_pruning_mpi
		cat test3.txt | mpirun -n 6 ./knapsack_dp_mpi
		cat test3.txt | mpirun -n 6 ./knapsack_pruning_mpi
comr:
		mpicc knapsack_dp_mpi.c -o knapsack_dp_mpi
		mpicc knapsack_pruning_mpi.c -o knapsack_pruning_mpi
		./generator 100 20 > ramdomtest.txt
		cat ramdomtest.txt | mpirun -n 6 ./knapsack_dp_mpi
		cat ramdomtest.txt | mpirun -n 6 ./knapsack_pruning_mpi
		rm ramdomtest.txt