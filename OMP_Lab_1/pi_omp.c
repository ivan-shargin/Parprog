// set number of threads through enviroment variable OMP_NUM_THREADS before running programm 
#include <stdio.h>
#include <omp.h>
static long num_steps = 100000000;
double step;
int main ()
{
	  double pi, sum = 0.0;
	  double start_time, run_time;
      int nthreads = 0;

	  step = 1.0/(double) num_steps;

        	 
	  start_time = omp_get_wtime();
      #pragma omp parallel
      {
        if (omp_get_thread_num() == 0){
            nthreads = omp_get_num_threads();
        } 
        #pragma omp for reduction(+:sum)
        for (int i = 1;i <= num_steps; i++){
            double x = (i-0.5)*step;
            sum = sum + 4.0/(1.0+x*x);
        } 
      }

	  pi = step * sum;
	  run_time = omp_get_wtime() - start_time;
	  printf("\n pi using %ld steps is %lf in %lf seconds\n ",num_steps,pi,run_time);
      printf("actual number of threads was %d\n", nthreads);
}	  