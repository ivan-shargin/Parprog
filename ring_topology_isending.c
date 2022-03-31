#include <mpi.h>
#include <stdio.h>
int main(int argc, char *argv[])
{
int size, rank, next, prev, buf=177, rbuf, tag = 777;
MPI_Request reqs[2];
MPI_Status stats[2];
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD, &size);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
prev = rank - 1;
if (prev < 0)
    prev = size - 1;
next = rank + 1;
if (next > size - 1)
    next = 0;

MPI_Irecv(&rbuf, 1, MPI_INT, prev, tag, MPI_COMM_WORLD, &reqs[1]);
printf("rank %d receive %d from %d  \n", rank, rbuf, prev);
MPI_Isend(&buf, 1, MPI_INT, next, tag, MPI_COMM_WORLD, &reqs[0]);
printf("rank %d send %d to %d  \n", rank, buf, next);

MPI_Waitall(2, reqs, stats);
MPI_Finalize();
}
