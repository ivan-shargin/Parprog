#include <stdio.h>
#include <mpi.h>
int main(int argc, char ** argv)
{
    int size, rank, tag = 77;
    double recv_buf;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
        double t1 = MPI_Wtime();
        MPI_Send(&t1, 1, MPI_DOUBLE, 1, tag, MPI_COMM_WORLD);
    }
    if (rank == 1)
    {
        MPI_Recv(&recv_buf, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, 0);
        double t2 = MPI_Wtime();
        double t = t2 - recv_buf;
        printf("time for lag is %f\n", t);
    }
    MPI_Finalize();
}
