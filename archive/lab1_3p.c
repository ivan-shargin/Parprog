#include <stdio.h>
#include <math.h>
#include <stdlib.h>

const int M = 1000;
const int K = 100;
const int len = (M + 1) * (K + 1);
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
    for (int j = 0;j < len; j++)
    {
        // f[j] = sin(2 * 3.14 * j / M) / 100;
        f[j] = 0;
    }

    for (int m = 0;m <= M; m++)
    {
        // U0[m] = (m * m - M * m + M * M / 4) / 1000.0;
        U0[m] = 0;
        Solution[m] = U0[m];
    }

    for (int m = 0;m <= 100; m++)
    {
        // U0[m] = (m * m - M * m + M * M / 4) / 1000.0;
        U0[m] = sin(2 * 3.14 * m / 100);
        Solution[m] = U0[m];
    }

    for (int k = 0;k < K; k++)
        for (int m = 0; m < M + 1; m++)
        {
            if ((m != 0) && (m != M))
                {
                    int i1 = ind(k, m - 1, K, M);
                    int i2 = ind(k, m + 1, K, M);
                    int i3 = ind(k + 1, m, K, M);
                    int i0 = ind(k, m, K, M);
                    Solution[i3] = (Solution[i1] + Solution[i2]) * 0.5 +
                    C * 0.5 * (Solution[i1] - Solution[i2]) + tau * f[i0];
                }
            if (m == 0)
                {
                    int i1 = ind(k, m, K, M);
                    int i2 = ind(k, m + 2, K, M);
                    int i3 = ind(k + 1, m, K, M);
                    Solution[i3] = Solution[i1] +
                    C * 0.5 * (Solution[i1] - Solution[i2]) + tau * f[i1];
                }
            if (m == M)
                {
                    int i1 = ind(k, m - 2, K, M);
                    int i2 = ind(k, m, K, M);
                    int i3 = ind(k + 1, m, K, M);
                    Solution[i3] = Solution[i2] +
                    C * 0.5 * (Solution[i1] - Solution[i2]) + tau * f[i2];
                }
        }

    FILE *file = NULL;
    file = fopen("solution.txt", "wb");
    if (file == NULL)
    {
        printf("Can't open the output file!");
        getchar();
        return -1;
    }

    for (int j = 0; j < len; j++)
    {
        fprintf(file, "%Lf", Solution[j]);
        fprintf(file, "%c", ',');
    }
    fclose(file);

    return 0;
}

int ind(int k, int m, int K, int M)
{
    return k * (M + 1) + m;
}
