#include <stdio.h>
#include <mpi.h>
int main ( int argc , char** argv )
{
int ntasks , mytask ;
MPI_Init (&argc , &argv ) ;
MPI_Comm_size ( MPI_COMM_WORLD , &ntasks ) ;
MPI_Comm_rank ( MPI_COMM_WORLD , &mytask ) ;
printf ( "Hello world from task %d of %d \n" , mytask , ntasks ) ;
MPI_Finalize ( ) ;
return 0 ;
}