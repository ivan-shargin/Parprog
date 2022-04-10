#include <mpi.h>
int ind(int k, int m, int K, int M)
{
    return k * (M + 1) + m;
}

int template(double* f, double* Solution, double C, double tau, int K, int M, int k, int m)
{
    int i1 = ind(k, m - 1, K, M);
    int i2 = ind(k, m + 1, K, M);
    int i3 = ind(k + 1, m, K, M);
    int i0 = ind(k, m, K, M);
    Solution[i3] = (Solution[i1] + Solution[i2]) * 0.5 +
    C * 0.5 * (Solution[i1] - Solution[i2]) + tau * f[i0];
}

int template_left(double* f, double* Solution, double C, double tau, int K, int M, int k, int m)
{
    int i1 = ind(k, m, K, M);
    int i2 = ind(k, m + 2, K, M);
    int i3 = ind(k + 1, m, K, M);
    Solution[i3] = Solution[i1] +
    C * 0.5 * (Solution[i1] - Solution[i2]) + tau * f[i1];
}

int template_right(double* f, double* Solution, double C, double tau, int K, int M, int k, int m)
{
    int i1 = ind(k, m - 2, K, M);
    int i2 = ind(k, m, K, M);
    int i3 = ind(k + 1, m, K, M);
    Solution[i3] = Solution[i2] +
    C * 0.5 * (Solution[i1] - Solution[i2]) + tau * f[i2];
}

int computing_cycle(int size, int rank, double* f, double* Solution, double C, double tau, int K, int M)
{
    int j = 0, m = 0, k = 0, min_m = 0, max_m = M;
    // MPI_Request *send_req = malloc(sizeof(MPI_Request) * size);
    // MPI_Request *recv_req = malloc(sizeof(MPI_Request) * (size-1));
    MPI_Request send_req;
    MPI_Request recv_req[10];
    void *buf = malloc(sizeof(int) * (size-1));

    if (size > 1){
        max_m = (M + 1)/ (size - 1) * (rank + 1);
        min_m = (M + 1)/ (size - 1) * rank;
        if (max_m > M) max_m = M;
    }
    else{
        max_m = M;
        min_m = 0;
    }

    for (k = 0;k < K; k++){
        for (m = min_m; m <= max_m; m++){
            if (((m != min_m) && (m != max_m)) ||
            ((k == 0) && (m != 0) &&(m != M))){
                template(f, Solution, C, tau, K, M, k, m);
            }

            if (m == 0){
                template_left(f, Solution, C, tau, K, M, k, m);
            }

            if (m == M){
                template_right(f, Solution, C, tau, K, M, k, m);
            }
        }
        if (k == 0){
            return 0;
        }

         if (rank != 0)
        {
            MPI_Isend(&rank, 1, MPI_INT, 0, 7, MPI_COMM_WORLD,
            &send_req);
            MPI_Wait(&send_req, MPI_STATUSES_IGNORE);
        }
        else
        {
            for (j=1; j < size; j++)
            {
                MPI_Irecv(buf + (j-1) * sizeof(int), 1, MPI_INT, j, 7, MPI_COMM_WORLD,
                &recv_req[j-1]);
                MPI_Waitall(size - 1, recv_req, MPI_STATUSES_IGNORE);
            }
        }

        free(buf);
    }
}
