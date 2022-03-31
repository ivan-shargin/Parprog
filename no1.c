#include <stdio.h>
#include <mpi.h>
int main ( int argc , char** argv )
{
int commsize, myrank ;
MPI_Init (&argc , &argv ) ;
MPI_Comm_size ( MPI_COMM_WORLD , &myrank ) ;
MPI_Comm_rank ( MPI_COMM_WORLD , &commsize ) ;
printf ( "Communicator size=%d My rank=%d \n" , commsize , myrank ) ;
MPI_Finalize ( ) ;
return 0 ;
}
