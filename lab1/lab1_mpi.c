#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include "computing_cycle.h"

const int M = 1000;
const int K = 500;
const int len = 1001 * 501;
int ind(int k, int m, int K, int M);

int main (int argc,char **argv)
{
    double x0 = 0, xM_1 = 1;
    double U0[M + 1];
    double a = 0.05;
    double h = (xM_1 - x0) / M;
    double T = 1;
    double tau = T / K;
    double C = a * tau / h;
    double Solution[len];
    double Times[100];
    // double *Times;
    double f[len];
    int j = 0, m = 0, k = 0;

    FILE *time_file = NULL;

    for (j = 0;j < len; j++)
    {
        // f[j] = sin(2 * 3.14 * j / M) / 100;
        f[j] = 0;
    }

    for (m = 0;m <= M; m++)
    {
        // U0[m] = (m * m - M * m + M * M / 4) / 1000.0;
        U0[m] = 0;
        Solution[m] = U0[m];
    }

    for (m = 0;m <= 100; m++)
    {
        // U0[m] = (m * m - M * m + M * M / 4) / 1000.0;
        U0[m] = -(m - 100) * (m -100) / 100 + 100;
        Solution[m] = U0[m];
    }

    int size = 0;
    int rank = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
        time_file = fopen("time.txt", "wb");
        if (time_file == NULL)
        {
            printf("Can't open the output file!");
            getchar();
            return -1;
        }
        fprintf(time_file, "%d", size);
        fprintf(time_file, "%c", ',');
        fclose(time_file);
    }

    double t1 = MPI_Wtime();
    computing_cycle(size, rank, f, Solution, C, tau, K, M);
    double t2 = MPI_Wtime();

    Times[rank] = t2 - t1;
    printf("For rank = %d computing took %f\n", rank, Times[rank]);

    for(j = 0;j < size; j++)
    {
        if (rank == j)
        {
            time_file = fopen("time.txt", "ab");
            if (time_file == NULL)
            {
                printf("Can't open the output file!");
                getchar();
                return -1;
            }
            fprintf(time_file, "%f", Times[j]);
            fprintf(time_file, "%c", ',');
            fclose(time_file);
        }
    }

    MPI_Finalize();

    FILE *file = NULL;
    file = fopen("solution.txt", "wb");
    if (file == NULL)
    {
        printf("Can't open the output file!");
        getchar();
        return -1;
    }

    for (k = 0; k < K; k += 10)
    {
        for (m = 0; m <= M; m ++)
        {
            j = ind(k, m, K, M);
            fprintf(file, "%f", Solution[j]);
            fprintf(file, "%c", ',');
        }
    }
    fclose(file);

    return 0;
}


int ind(int k, int m, int K, int M)
{
    return k * (M + 1) + m;
}
