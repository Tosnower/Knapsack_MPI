all:
	mpicc knap_MPI.c -o Knap_MPI
	cat test30.txt | mpirun -n 4 ./Knap_MPi
de:
	mpicc knap_MPI.c -o Knap_MPI
	./generator 100 50 | mpirun -n 6 ./Knap_MPi
debug:
	mpicc -DDEBUG knap_MPI.c -o Knap_MPI
tmp:
		mpicc tmp.c -o tmp
		./generator 100 10 | mpirun -n 2 ./tmp
