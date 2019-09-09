#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <sys/time.h>

 int knapSack( long int C,  long int w[],  long int v[], int n);


uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}

int main(int argc, char *argv[]) {
    long int C;    /* capacity of backpack */
    int n;    /* number of items */
    int i;    /* loop counter */

    MPI_Init (&argc, &argv);
    int rank;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        scanf ("%d", &C);
        scanf ("%d", &n);
    }

    long int v[n], w[n];        /* value, weight */

    if (rank == 0) {
        for (i = 0; i < n; i++) {
            scanf ("%d %d", &v[i], &w[i]);
        }
    }


    uint64_t start = GetTimeStamp ();
    int ks = knapSack(C, w, v, n);

    if (rank == 0) {
        printf ("knapsack occupancy %ld\n", ks);
        printf ("Time: %ld us\n", (uint64_t) (GetTimeStamp() - start));
    }

    MPI_Finalize ();

    return 0;
}


/* PLACE YOUR CHANGES BELOW HERE */
//subblock size
#define ROW 1
#define COL 512

int min(int i, int j) {
    return (i<j) ? i : j;
}


long int solver(int n, long int c, int rows, int weight[rows], int profit[rows],
            int start, int rank, int size) {
    int recv_rank = (rank-1)%size;   //rank to receive data
    int send_rank = (rank+1)%size; // rank to send data
    if( start == 0 ) {      // deal with first block, since it doesn't receive data from any nodes

        long int total[rows][c];
        int i, j;

        for (j = 0; j < c; j += COL) {
            int cols = min(COL, c-j);
            int k;
            for (k = j; k < j + cols; k++) {
                if (weight[0] > k) {
                    total[0][k] = 0;
                } else {
                    total[0][k] = profit[0];
                }
            }
            //compute subblock
            for (i = 1; i < rows; i++) {
                for (k = j; k < j + cols; k++) {
                    //int ni = i+start;
                    if ( (k<weight[i]) ||
                        (total[i-1][k] >= total[i-1][k-weight[i]] + profit[i])) {
                        total[i][k] = total[i-1][k];
                    } else {
                        total[i][k] = total[i-1][k-weight[i]] + profit[i];
                    }
                }
            }
            //send last row to next node
            MPI_Send(&total[rows-1][j], cols, MPI_LONG, send_rank, j, MPI_COMM_WORLD);
        }

    } else {
//        printf("I run1");
        long int total[rows+1][c];  // use the first row to store the data from last node
        int i, j;
        for (j = 0; j < c; j += COL) {
            int cols = min(COL, c-j);
            int k;
            // receive data from last node
            MPI_Recv(&total[0][j], cols, MPI_LONG, recv_rank, j, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (i = 1; i <= rows; i++) {
                for (k = j; k < j + cols; k++) {
                    int ni = i-1;  //ni is index for weight and profit, notice the first row is the data from last node
                    if ( (k<weight[ni]) ||
                        (total[i-1][k] >= total[i-1][k-weight[ni]] + profit[ni])) {
                        total[i][k] = total[i-1][k];
                    } else {
                        total[i][k] = total[i-1][k-weight[ni]] + profit[ni];
                    }
                }
            }


            if (start + rows == n && j + cols == c) {
                //computer the last subblock of last ROW, print the final result
                printf("max profit: %d \n", total[rows][c-1]);
                return total[rows][c-1];
            } else if (start + rows != n){
                // if it is not last ROW, we need send data to next node.
                MPI_Send(&total[rows][j], cols, MPI_LONG, send_rank, j, MPI_COMM_WORLD);
            }
        }
    }
}



 int knapSack(long int C,  long int w[], long int v[], int n) {
    int size, rank,i;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        for (int k = 1; k < size; k++) {
            MPI_Send(&n, 1, MPI_INT, k, 0, MPI_COMM_WORLD);
            MPI_Send(&C, 1, MPI_LONG, k, 1, MPI_COMM_WORLD);
            MPI_Send(v, n, MPI_LONG, k, 2, MPI_COMM_WORLD);
            MPI_Send(w, n, MPI_LONG, k, 3, MPI_COMM_WORLD);
        }

    } else {
        MPI_Recv(&n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&C, 1, MPI_LONG, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        v = (int*) malloc(sizeof(long int) * n);
        w = (int*) malloc(sizeof(long int) * n);
        MPI_Recv(v, n, MPI_LONG, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(w, n, MPI_LONG, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    for (i = 0; i < n; i += ROW) {
        int rows  = min(ROW, n-i);
        int *weights = malloc(rows * sizeof(int));   //initial weights locally
        int *profits = malloc(rows * sizeof(int));   //initial weights locally
        for (int j = 0; j<rows; j++){
            weights[j] = w[i+j];
            profits[j] = v[i+j];
        }
        if ((i/ROW) % size == rank)
            solver(n, C, rows, weights, profits, i, rank, size);  //solve subblock
        free(weights);
        free(profits);
    }
}
