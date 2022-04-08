int computing_cycle(int size, int rank, double* f, double* Solution, double C, double tau, int K, int M)
{
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
}
