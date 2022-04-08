#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>

const int M = 1000;
const int K = 100;
const int len = 1001 * 101;
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
    double f[len];
    double arrTimes[100];

    FILE *Times = NULL;
    Times = fopen("Times.txt", "wb");
    if (Times == NULL)
    {
        printf("Can't open the output file!");
        getchar();
        return -1;
    }

    int size = 0;
    int rank = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double t1 = MPI_Wtime();
    int j = 0, m = 0, k = 0, min_m = 0, max_m = M;
    if (size > 1)
    {
        max_m = (M + 1)/ (size - 1) * (rank + 1);
        min_m = (M + 1)/ (size - 1) * rank;
        if (max_m > M) max_m = M;
    }
    else
    {
        max_m = M;
        min_m = 0;
    }

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

    for (k = 0;k < K; k++)
        for (m = min_m; m <= max_m; m++)
        {
            if ((m != min_m) && (m != max_m))
                {
                    int i1 = ind(k, m - 1, K, M);
                    int i2 = ind(k, m + 1, K, M);
                    int i3 = ind(k + 1, m, K, M);
                    int i0 = ind(k, m, K, M);
                    Solution[i3] = (Solution[i1] + Solution[i2]) * 0.5 +
                    C * 0.5 * (Solution[i1] - Solution[i2]) + tau * f[i0];
                }
            if (m == min_m)
                {
                    int i1 = ind(k, m, K, M);
                    int i2 = ind(k, m + 2, K, M);
                    int i3 = ind(k + 1, m, K, M);
                    Solution[i3] = Solution[i1] +
                    C * 0.5 * (Solution[i1] - Solution[i2]) + tau * f[i1];
                }
            if (m == max_m)
                {
                    int i1 = ind(k, m - 2, K, M);
                    int i2 = ind(k, m, K, M);
                    int i3 = ind(k + 1, m, K, M);
                    Solution[i3] = Solution[i2] +
                    C * 0.5 * (Solution[i1] - Solution[i2]) + tau * f[i2];
                }
        }
    double t2 = MPI_Wtime();
    printf("For rank = %d computing took %f\n", rank, t2-t1);
    fprintf(Times, "%f", t2-t1);
    fprintf(Times, "%c", ',');
    MPI_Finalize();
    // double t2 = MPI_Wtime();
    // printf("For rank = %d computing took %f\n", rank, t2-t1);

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
    fclose(Times);

    return 0;
}

int ind(int k, int m, int K, int M)
{
    return k * (M + 1) + m;
}
