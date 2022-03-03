#include <stdio.h>
#include <mpi.h>
int main ( int argc , char** argv )
{

	int N = 1000000; 
	double sum, sum_send, sum_recv; 
	double t1, t2; //time
	
	int i , rank , n_rank , count , start , stop;
    
	MPI_Init(&argc , &argv ) ;
	MPI_Comm_rank ( MPI_COMM_WORLD , &rank ) ;
	MPI_Comm_size ( MPI_COMM_WORLD , &n_rank ) ;
	count = N / n_rank ;
    
	start = (rank + 1) * count ;
    
	stop = start + count ;
	sum = 0;
	sum_send = 0;
	
	
	if ( rank != 0 ) 
	{
		t1 = MPI_Wtime();
		for(i = start; i < stop; i++)
		{
            sum_send  +=  (double) 1/i;
		}
		
		MPI_Send (&sum_send , 1 , MPI_DOUBLE , 0 , 0 , MPI_COMM_WORLD );
		t2 = MPI_Wtime();
		printf("rank %d time mpi: %f\n",rank,  t2- t1);
	} else 
	{
		t1 = MPI_Wtime();
		for(i = start; i < stop; i++)
		{
            sum  +=  (double) 1/i;
		}
		for ( i =1; i < n_rank; i++ ) 
		{
			MPI_Recv(&sum_recv , 1 , MPI_DOUBLE, i , 0 , MPI_COMM_WORLD, 0 ) ;
			sum += sum_recv;
		}
		t2 = MPI_Wtime();
		printf("rank %d time mpi: %f\n",rank,  t2- t1);
		
	}
	if(rank == 0)
		printf ( "sum =  %f \n" , sum) ;

	MPI_Finalize ( ) ;
	
	return 0;
}