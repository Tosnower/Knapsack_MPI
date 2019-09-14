pr:
	mpicc knapsack_pruning_mpi.c -o knapsack_pruning_mpi
	./generator 10000000 35 | mpirun -n 6 ./knapsack_pruning_mpi
prc:
	mpicc knapsack_pruning_mpi.c -o knapsack_pruning_mpi
	cat test300100.txt | mpirun -n 6 ./knapsack_pruning_mpi

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
		cat test30.txt | mpirun -n 6 ./knapsack_dp_mpi
		cat test30.txt | mpirun -n 6 ./knapsack_pruning_mpi
comr:
		mpicc knapsack_dp_mpi.c -o knapsack_dp_mpi
		mpicc knapsack_pruning_mpi.c -o knapsack_pruning_mpi
		./generator 100 20 > ramdomtest.txt
		cat ramdomtest.txt | mpirun -n 6 ./knapsack_dp_mpi
		cat ramdomtest.txt | mpirun -n 6 ./knapsack_pruning_mpi
		rm ramdomtest.txt
final:
		mpicc tosnower-knapsack_finalv2.c -o knapsack_finalv
		./generator 5000 500 | mpirun -n 6 ./knapsack_finalv
seq:
		mpicc knap_sequential.c -o knap_sequentialv
		mpicc knapsack_dp_mpi.c -o knapsack_dp_mpi
		mpicc knapsack_pruning_mpi.c -o knapsack_pruning_mpi
		./generator 10000 100 > test10020.txt
		./generator 40000 100 > test10030.txt
		./generator 80000 100 > test10040.txt
		./generator 160000 100 > test10050.txt
		./generator 320000 100 > test10060.txt
		cat test10020.txt | mpirun -n 6 ./knapsack_pruning_mpi
		cat test10030.txt | mpirun -n 6 ./knapsack_pruning_mpi
		cat test10040.txt | mpirun -n 6 ./knapsack_pruning_mpi
		cat test10050.txt | mpirun -n 6 ./knapsack_pruning_mpi
		cat test10060.txt | mpirun -n 6 ./knapsack_pruning_mpi
		cat test10020.txt | mpirun -n 6 ./knapsack_dp_mpi
		cat test10030.txt | mpirun -n 6 ./knapsack_dp_mpi
		cat test10040.txt | mpirun -n 6 ./knapsack_dp_mpi
		cat test10050.txt | mpirun -n 6 ./knapsack_dp_mpi
		cat test10060.txt | mpirun -n 6 ./knapsack_dp_mpi
		rm test10020.txt
		rm test10030.txt
		rm test10040.txt
		rm test10050.txt
		rm test10060.txt