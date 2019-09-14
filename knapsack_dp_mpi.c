#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <sys/time.h>

long int knapSack( long int C,  long int w[],  long int v[], int n);


uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}

int main(int argc, char *argv[]) {
    long int C;    /* capacity of backpack */
    int n;    /* number of items */
    int i;    /* loop counter */

    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        scanf("%ld", &C);
        scanf("%d", &n);
    }

    MPI_Bcast(&C, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    long int v[n], w[n];        /* value, weight */

    if (rank == 0) {
        for (i = 0; i < n; i++) {
            scanf("%ld %ld", &v[i], &w[i]);
        }

    }
    MPI_Bcast(v, n, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(w, n, MPI_LONG, 0, MPI_COMM_WORLD);


    uint64_t start = GetTimeStamp();
    long int ks = knapSack(C, w, v, n);

    if (rank == 0) {
        printf("knapsack occupancy %ld\n", ks);
        printf("Time: %ld us\n", (uint64_t) (GetTimeStamp() - start));
    }

    MPI_Finalize();

    return 0;
}


/* PLACE YOUR CHANGES BELOW HERE */
//subblock size


long int knapSack(long int C, long int w[], long int v[], int n) {
    ++C;
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rowSize = n/size+1;
    int colSize = C/size+1;
    int i = 0;
    if (i==0 && rank==0) {
        int blockRow = (rowSize < n-i)?rowSize:n-i;
        long int blockWeight[blockRow];
        long int blockValue[blockRow];
        for (int j = 0; j < blockRow; j++) {
            blockWeight[j] = w[i + j];
            blockValue[j] = v[i + j];
        }
        long int tmp[blockRow][C];
        for (int j = 0; j < C; j += colSize) {
            int blockCol = (colSize<C-j)? colSize:C-j;
            for (int k = j; k < j + blockCol; k++) {
                if (blockWeight[0] > k) {
                    tmp[0][k] = 0;
                } else {
                    tmp[0][k] = blockValue[0];
                }
            }
            for (int x = 1; x < blockRow; x++) {
                for (int k = j; k < j + blockCol; k++) {
                    if (k < blockWeight[x]) {
                        tmp[x][k] = tmp[x - 1][k];
                    } else if ( tmp[x - 1][k - blockWeight[x]] + blockValue[x] > tmp[x - 1][k]){
                        tmp[x][k] = tmp[x - 1][k - blockWeight[x]] + blockValue[x];
                    } else{
                        tmp[x][k] = tmp[x - 1][k];
                    }
                }
            }
            MPI_Send(&tmp[blockRow - 1][j], blockCol, MPI_LONG, (rank + 1) % size, j, MPI_COMM_WORLD);
        }
        i+=rowSize;
    }
    for (; i < n; i += rowSize) {
        int blockRow = (rowSize < n-i)?rowSize:n-i;
        long int blockWeight[blockRow];
        long int blockValue[blockRow];
        for (int j = 0; j < blockRow; j++) {
            blockWeight[j] = w[i + j];
            blockValue[j] = v[i + j];
        }
        if ((i / rowSize) % size == rank) {
            long int tmp[blockRow + 1][C];
            for (int j = 0; j < C; j += colSize) {
                int blockCol = (colSize<C-j)? colSize:C-j;
//                int blockCol = min(colSize, C - j);
                MPI_Recv(&tmp[0][j], blockCol, MPI_LONG, (rank - 1) % size, j, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int x = 1; x <= blockRow; x++) {
                    for (int k = j; k < j + blockCol; k++) {
                        if (k < blockWeight[x-1]) {
                            tmp[x][k] = tmp[x - 1][k];
                        } else if ( tmp[x - 1][k - blockWeight[x-1]] + blockValue[x-1] > tmp[x - 1][k]){
                            tmp[x][k] = tmp[x - 1][k - blockWeight[x-1]] + blockValue[x-1];
                        } else{
                            tmp[x][k] = tmp[x - 1][k];
                        }
                    }
                }
                if (i + blockRow == n && j + blockCol == C) {
                    MPI_Send(&tmp[blockRow][C - 1], 1, MPI_LONG, 0, 0, MPI_COMM_WORLD);
                } else if (i + blockRow != n) {
                    MPI_Send(&tmp[blockRow][j], blockCol, MPI_LONG, (rank + 1) % size, j, MPI_COMM_WORLD);
                }
            }
        }
    }
    long res = 0;
    if (rank == 0) {
        MPI_Status status;
        int tag;
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        tag = status.MPI_TAG;
        int src = status.MPI_SOURCE;
        if (tag == 0) {
            MPI_Recv(&res, n + 1, MPI_LONG, src, tag, MPI_COMM_WORLD, &status);
            return res;
        }
    }
    return 0;
}