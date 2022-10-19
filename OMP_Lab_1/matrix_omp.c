#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "print_matrix.h"
#include "transpose.h"
#include "convolute.h"
#include <omp.h>

int main()
{
    const int N = 500;

    // double *A = (double *)malloc(sizeof(double) * N * N);
    // double *B = (double *)malloc(sizeof(double) * N * N);
    double A[N*N];
    double B[N*N];

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            A[i * N + j] = i * N + j;
            if (i == j){
                B[i * N + j] = 10.0;
            } else {
                B[i * N + j] = 0.0;
            }
        }
    }

    // double *C = (double *)malloc(sizeof(double) * N * N);
    double C[N*N];
    double *B_tr = (double *)malloc(sizeof(double) * N * N);
    // free(B);
    transpose(B, B_tr, N);
    #pragma omp parallel firstprivate(A, B_tr, N)
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            C[i * N + j] = convolute(A + i * N, B_tr + j * N, N);
        }
    }

    free(B_tr);
    // free(A);
    // free(C);


    
    // double A_test[4] = {10.0, 0.0, 0.0, 10.0};
    // double B_test[4] = {1.0, 2.0, 3.0, 4.0};
    // int n = 2;
    // print_matrix(A_test, n);
    // print_matrix(B_test, n);

    // double C_test[4];
    // double *B_test_tr = (double *)malloc(sizeof(double) * n * n);
    // transpose(B_test, B_test_tr, n);
    // #pragma omp parallel firstprivate(A_test, B_test_tr, n)
    // for (int i = 0; i < n; i++){
    //     for (int j = 0; j < n; j++){
    //         C_test[i * n + j] = convolute(A_test + i * n, B_test_tr + j * n, n);
    //     }
    // }

    // print_matrix(C_test, n);
    // free(B_test_tr);

    return 0;    
}