#include <stdio.h>
#include <mpi.h>
int main ( int argc , char** argv )
{

	int N = 1000000000; 
	double e, e_send, e_recv; 
	double t1, t2; //time
	
	int i , rank , n_rank , count , start , stop;
    
	MPI_Init(&argc , &argv ) ;
	MPI_Comm_rank ( MPI_COMM_WORLD , &rank ) ;
	MPI_Comm_size ( MPI_COMM_WORLD , &n_rank ) ;
	count = N / n_rank ;
    
	start = rank * count ;
    
	stop = start + count ;
	e = 1;
	e_send = 1;
	
	
	if ( rank != 0 ) 
	{
		t1 = MPI_Wtime();
		for(i = start; i < stop; i++)
		{
            e_send  *=  (1 + (double) 1/N);
		}
		
		MPI_Send (&e_send , 1 , MPI_DOUBLE , 0 , 0 , MPI_COMM_WORLD );
		t2 = MPI_Wtime();
		printf("rank %d time mpi: %f\n",rank,  t2- t1);
	} else 
	{
		t1 = MPI_Wtime();
		for(i = start; i < stop; i++)
		{
            e  *= (1+ (double) 1/N);
		}
		for ( i =1; i < n_rank; i++ ) 
		{
			MPI_Recv(&e_recv , 1 , MPI_DOUBLE, i , 0 , MPI_COMM_WORLD, 0 ) ;
			e *= e_recv;
		}
		t2 = MPI_Wtime();
		printf("rank %d time mpi: %f\n",rank,  t2- t1);
		
	}
	if(rank == 0)
		printf ( "e =  %f \n" , e) ;

	MPI_Finalize ( ) ;
	
	return 0;
}