#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "print_matrix.h"
#include "transpose.h"
#include "convolute.h"
#include <omp.h>

int main()
{
    const int k = 5;
    int N = pow(2, k);

    double *A = (double *)calloc(sizeof(double), N * N);
    double *B = (double *)calloc(sizeof(double), N * N);
  
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            A[i * N + j] = (i * N + j) / N / N;
            if (i == j){
                B[i * N + j] = 10.0;
            } else {
                B[i * N + j] = 0.0;
            }
        }
    }

    double *C = (double *)calloc(sizeof(double), N * N);  

    //direct strassen
      
    

    return 0;    
}