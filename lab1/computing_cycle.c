#include <mpi.h>
#include <stdio.h>
int ind(int k, int m, int width)
{
    return k * width + m;
}

int template(double *Solution, double *f, int i0, int i1, int i2, int i3, double C, double tau)
{
    Solution[i3] = (Solution[i1] + Solution[i2]) * 0.5 +
    C * 0.5 * (Solution[i1] - Solution[i2]) + tau * f[i0];
}

int template_left(double *Solution, double *f, int i1, int i2, int i3, double C, double tau)
{
    Solution[i3] = Solution[i1] +
    C * 0.5 * (Solution[i1] - Solution[i2]) + tau * f[i1];
}

int template_right(double *Solution, double *f, int i1, int i2, int i3, double C, double tau)
{
    Solution[i3] = Solution[i2] +
    C * 0.5 * (Solution[i1] - Solution[i2]) + tau * f[i2];
}


int computing_cycle(int size, int rank, double* f, double* U0, double* Solution, double C, double tau, int K, int M)
{
    int j = 0, m = 0, k = 0, min_m = 0, max_m = M, width = M;
    MPI_Request send_reqs[2];
    MPI_Request recv_reqs[2];
    int tag = 7;
    int prev = rank - 1;
    int next = rank + 1;

    if (size > 1){
        max_m = M / size * (rank + 1);
        min_m = M / size * rank;
        if ((rank == size - 1) && (max_m < M)) max_m = M;
        width = max_m - min_m;
    }
    else{
        max_m = M;
        min_m = 0;
        width = M;
    }
    printf("process with rank = %d has width = %d min_m = %d max_m = %d\n", rank, width, min_m, max_m);

    for (m = min_m; m < max_m; m++){
        int m1 = m - min_m;
        Solution[m1] = U0[m1];
    }

    for (m = min_m; m < max_m; m++){
        if ((m != 0) && (m != M - 1)){
            int m1 = m - min_m;
            int i1 = ind(0, m1 - 1, width);
            int i2 = ind(0, m1 + 1, width);
            int i3 = ind(1, m1, width);
            int i0 = ind(0, m1, width);
            Solution[i3] = (U0[m1 - 1] + U0[m1 + 1])* 0.5 + C * 0.5 * (U0[m1 - 1] - U0[m1 + 1]) + tau * f[i0];
        }

        if (m == 0){
            int i1 = ind(0, 0, width);
            int i2 = ind(0, 2, width);
            int i3 = ind(1, 0, width);
            Solution[i3] = U0[i1] +
            C * 0.5 * (U0[i1] - U0[i2]) + tau * f[i1];
        }

        if (m == M - 1){
            int m1 = m - min_m;
            int i1 = ind(0, m1 - 2, width);
            int i2 = ind(0, m1, width);
            int i3 = ind(1, m1, width);
            Solution[i3] = U0[i2] +
            C * 0.5 * (U0[i1] - U0[i2]) + tau * f[i2];
        }
    }

    if (prev >= 0){
        MPI_Isend(Solution, 1, MPI_DOUBLE,
        prev, tag, MPI_COMM_WORLD, send_reqs);
    }

    if (next < size){
        MPI_Isend(Solution + width - 1, 1, MPI_DOUBLE,
        next, tag, MPI_COMM_WORLD, send_reqs + 1);
    }

    for (k = 1;k < K - 1; k++){

        double left, right;
        if (prev >= 0){
            MPI_Irecv(&left, 1, MPI_DOUBLE, prev, tag,
            MPI_COMM_WORLD, recv_reqs);
        }

        if (next < size){
            MPI_Irecv(&right, 1, MPI_DOUBLE, next, tag,
            MPI_COMM_WORLD, recv_reqs + 1);
        }

        for (m = min_m; m < max_m; m++){
            int m1 = m - min_m;
            if ((m != min_m) && (m != max_m)){
                int i1 = ind(k, m1 - 1, width);
                int i2 = ind(k, m1 + 1, width);
                int i3 = ind(k + 1, m1, width);
                int i0 = ind(k, m1, width);
                Solution[i3] = (Solution[i1] + Solution[i2]) * 0.5 +
                C * 0.5 * (Solution[i1] - Solution[i2]) + tau * f[i0];
                // template(f, Solution, i0, i1, i2, i3, C, tau);
            }

            if (m == 0){
                int i1 = ind(k, m1, width);
                int i2 = ind(k, m1 + 2, width);
                int i3 = ind(k + 1, m1, width);
                Solution[i3] = Solution[i1] +
                C * 0.5 * (Solution[i1] - Solution[i2]) + tau * f[i1];
                // template_left(f, Solution, i1, i2, i3, C, tau);
            }

            if (m == M-1){
                int i1 = ind(k, m1 - 2, width);
                int i2 = ind(k, m1, width);
                int i3 = ind(k + 1, m1, width);
                Solution[i3] = Solution[i2] +
                C * 0.5 * (Solution[i1] - Solution[i2]) + tau * f[i2];
                // template_right(f, Solution, i1, i2, i3, C, tau);
            }
        }

        if (prev >= 0){
            MPI_Wait(recv_reqs, MPI_STATUSES_IGNORE);
            int i2 = ind(k, 1, width);
            int i3 = ind(k + 1, 0, width);
            int i0 = ind(k, 0, width);
            Solution[i3] = (left + Solution[i2]) * 0.5 +
            C * 0.5 * (left - Solution[i2]) + tau * f[i0];
            MPI_Isend(&Solution[i3], 1, MPI_DOUBLE, prev, tag,
            MPI_COMM_WORLD, send_reqs);
        }

        if (next < size){
            MPI_Wait(recv_reqs + 1, MPI_STATUSES_IGNORE);
            int i1 = ind(k, width - 2, width);
            int i3 = ind(k + 1, width - 1, width);
            int i0 = ind(k, width - 1, width);
            Solution[i3] = (Solution[i1] + right) * 0.5 +
            C * 0.5 * (Solution[i1] - right) + tau * f[i0];
            MPI_Isend(&Solution[i3], 1, MPI_DOUBLE, next, tag,
            MPI_COMM_WORLD, send_reqs + 1);
        }
    }
}
