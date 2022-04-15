#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include "computing_cycle.h"

const int M = 10000;
const int K = 1000;
const double x0 = 0, xM = 1;
const double a = 0.5;
const double T = 1;
const int D_k = 100;

int main (int argc,char **argv)
{
    double U0[M];
    double h = (xM - x0) / (M - 1);
    double tau = T / (K - 1);
    double C = a * tau / h;
    int j = 0, i = 0, m = 0, k = 0;
    int len = M * K;
    int min_m = 0, max_m = M, width = M;
    int size = 0, rank = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (size > 1){
        max_m = M / size * (rank + 1);
        min_m = M / size * rank;
        if (max_m > M) max_m = M;
        width = max_m - min_m;
    }
    else{
        max_m = M;
        min_m = 0;
        width = M;
    }

    for (m = 0;m < M; m++){
        U0[m] = 0;
    }

    for (m = 0;m <= 100; m++){
        U0[m] = -(m - 100) * (m -100) / 100;
    }

    double *f = (double *) malloc(sizeof(double) * width * K);
    for (j = 0;j < width * K; j++){
        f[j] = 0;
    }

    double *Solution = (double*) malloc(width * K * sizeof(double));
    double t1 = MPI_Wtime();
    computing_cycle(size, rank, f, &U0[min_m], Solution, C, tau, K, M);
    double t2 = MPI_Wtime();

    double time = t2 - t1;
    printf("For rank = %d computing took %f\n", rank, time);

    int root = 0;

    if (rank == root){
        double *output = (double*) malloc(M * K / D_k * sizeof(double));
        FILE *file = NULL;
        file = fopen("output.bin", "wb");
        if (file == NULL){
            printf("Can't open the output file!");
            getchar();
            return -1;
        }

        for(k = 0;k < K / D_k;k++){
            for(j = 1;j < size;j++){
                MPI_Recv(output + k * M + j * width, );
            }
        }

        fclose(file);
        free(output);
    }


    free(Solution);
    free(f);
    MPI_Finalize();
    return 0;
}
