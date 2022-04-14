#include <mpi.h>
int ind(int k, int m, int M)
{
    return k * M + m;
}

int template(double *Solution, double *f, int i0, int i1, int i2, int i3, double C, double tau)
{
    Solution[i3] = (Solution[i1] + Solution[i2]) * 0.5 +
    C * 0.5 * (Solution[i1] - Solution[i2]) + tau * f[i0];
}

int template_left(double* Solution, double *f, int i1, int i2, int i3, double C, double tau)
{
    Solution[i3] = Solution[i1] +
    C * 0.5 * (Solution[i1] - Solution[i2]) + tau * f[i1];
}

int template_right(double* Solution, double *f, int i1, int i2, int i3, double C, double tau)
{
    Solution[i3] = Solution[i2] +
    C * 0.5 * (Solution[i1] - Solution[i2]) + tau * f[i2];
}

int computing_cycle(int size, int rank, double* f, double* U0, double* Solution, double C, double tau, int K, int M)
{
    int j = 0, m = 0, k = 0, min_m = 0, max_m = M - 1, width = M;
    MPI_Request send_reqs[2];
    MPI_Request recv_reqs[2];

    if (size > 1){
        max_m = M / (size - 1) * (rank + 1);
        min_m = M / (size - 1) * rank;
        if (max_m >= M) max_m = M - 1;
        width = max_m - min_m + 1;
    }
    else{
        max_m = M - 1;
        min_m = 0;
        width = M;
    }

    for (m = min_m; m <= max_m; m++){
        int m1 = m - min_m;
        Solution[m1] = U0[m1];
    }

    for (m = min_m; m <= max_m; m++){
        int m1 = m - min_m;
        int i1 = ind(0, m1 - 1, width);
        int i2 = ind(0, m1 + 1, width);
        int i3 = ind(1, m1, width);
        int i0 = ind(0, m1, width);
        Solution[i3] = (U0[i1] + U0[i2]) * 0.5 +
        C * 0.5 * (U0[i1] - U0[i2]) + tau * f[i0];
    }
    int tag = 7;
    int prev = rank - 1;
    int next = rank + 1;

    if (prev >= 0){
    MPI_Isend(Solution, 1, MPI_DOUBLE,
    int prev, tag, MPI_COMM_WORLD, send_reqs);
    }

    if (next < size){
    MPI_Isend(Solution + width - 1, 1, MPI_DOUBLE,
    int next, tag, MPI_COMM_WORLD, send_reqs + 1);
    }

    for (k = 0;k < K - 1; k++){
        double left, right;
        if (prev >= 0){
            MPI_Irecv(&left, 1, MPI_DOUBLE, prev, tag,
            MPI_COMM_WORLD, recv_reqs);
        }

        if (next < size){
            MPI_Irecv(&right, 1, MPI_DOUBLE, next, tag,
            MPI_COMM_WORLD, recv_reqs + 1);
        }I

        for (m = min_m; m <= max_m; m++){
            int m1 = m - min_m;
            if ((m != min_m) && (m != max_m)){
                int i1 = ind(k, m1 - 1, width);
                int i2 = ind(k, m1 + 1, width);
                int i3 = ind(k + 1, m1, width);
                int i0 = ind(k, m1, width);
                template(f, Solution, i0, i1, i2, i3, C, tau);
            }

            if (m == 0){
                int i1 = ind(k, m1, width);
                int i2 = ind(k, m1 + 2, width);
                int i3 = ind(k + 1, m1, width);
                template_left(f, Solution, i1, i2, i3, C, tau);
            }

            if (m == M){
                int i1 = ind(k, m1 - 2, K, width);
                int i2 = ind(k, m1, width);
                int i3 = ind(k + 1, m1, width);
                template_right(f, Solution, i1, i2, i3, C, tau);
            }

            if ((m == min_m) && (m != 0)){

            }
        }
    }
}
