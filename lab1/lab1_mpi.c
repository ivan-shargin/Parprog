#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include "computing_cycle.h"

const int M = 1000;
const int K = 100;

int main (int argc,char **argv)
{
    double x0 = 0, xM_1 = 1;
    double U0[M + 1];
    double a = 0.5;
    double h = (xM_1 - x0) / M;
    double T = 1;
    double tau = T / K;
    double C = a * tau / h;
    int j = 0, m = 0, k = 0;
    int len = (M + 1) * (K + 1);

    for (j = 0;j < len; j++)
    {
        f[j] = 0;
    }

    for (m = 0;m <= M; m++)
    {
        U0[m] = 0;
    }

    for (m = 0;m <= 100; m++)
    {
        U0[m] = -(m - 100) * (m -100) / 100 + 100;
    }

    int size = 0;
    int rank = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (size > 1){
        int width = (M + 1)/ (size - 1) + 1
        max_m = (M + 1)/ (size - 1) * (rank + 1);
        min_m = (M + 1)/ (size - 1) * rank;
        if (max_m > M) max_m = M;
    }
    else{
        width = M + 1;
        max_m = M;
        min_m = 0;
    }

    double Solution = malloc(sizeof(double) * width);
    double t1 = MPI_Wtime();
    computing_cycle(size, rank, f, &U0[min_m], Solution, C, tau, K, M);
    double t2 = MPI_Wtime();

    Times[rank] = t2 - t1;
    printf("For rank = %d computing took %f\n", rank, Times[rank]);

    if (rank == 0){
        FILE *file = NULL;
        file = fopen("solution.txt", "wb");
        if (file == NULL){
            printf("Can't open the output file!");
            getchar();
            return -1;
        }

        for (k = 0; k < K; k += 10){
            for (m = 0; m <= M; m ++){
                j = ind(k, m, K, M);
                fprintf(file, "%f", Solution[j]);
                fprintf(file, "%c", ',');
            }
        }
        fclose(file);
        free(Solution);
    }
    MPI_Finalize();




    return 0;
}
