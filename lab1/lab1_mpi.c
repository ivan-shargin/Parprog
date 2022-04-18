#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include "computing_cycle.h"

const int M = 1000;
const int K = 1000;
const double x0 = 0, xM = 1;
const double a = 0.2;
const double T = 1;
const int D_k = 10;

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

    for (m = 3. * M / 10.; m < 7. * M / 10.; m++){
        U0[m] = sin(2 * 3.14 * m / 500);
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
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == root){
        double *output = (double*) malloc(M * K / D_k * sizeof(double));

        for(i = 0;i < K / D_k;i ++){
            k = i * D_k;
            MPI_Gather(Solution + k * width, width, MPI_DOUBLE,
            output + i * M, width, MPI_DOUBLE, root, MPI_COMM_WORLD);
        }

        FILE *file = NULL;
        file = fopen("output.bin", "wb");
        if (file == NULL){
            printf("Can't open the output file!");
            getchar();
            return -1;
        }

        for(i = 0;i < K / D_k * M; i++){
            fprintf(file, "%f", output[i]);
            fprintf(file, "%c", ',');
        }

        free(output);
        fclose(file);

    }else{
        for(i = 0;i < K / D_k;i ++){
            k = i * D_k;
            MPI_Gather(Solution + k * width, width, MPI_DOUBLE,
            NULL, width , MPI_DOUBLE, root, MPI_COMM_WORLD);
        }

    }



    free(Solution);
    free(f);
    MPI_Finalize();
    return 0;
}
