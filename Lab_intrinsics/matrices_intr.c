#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <emmintrin.h>

int print_matrix(double *Matrix, int n){
    printf("\n");
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            printf("%8.1lf", Matrix[i * n + j]);
        }
        printf("\n");
    }
    return 0;
}

int transpose(double *Matrix, double *Matrix_tr, int n){    
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            Matrix_tr[j * n + i] = Matrix[i * n + j];
        }        
    }
    return 0;
}

double convolute(double *A, double *B, int n){
    double sum = 0.0;
    __m128d _sum, _Aelem, _Belem, _mul;
    _sum = _mm_set_pd(0.0, 0.0);
    for (int i = 0; i < n; i+=2){
        _Aelem = _mm_set_pd(A[i], A[i+1]);
        _Belem = _mm_set_pd(B[i], B[i+1]);
        _mul = _mm_mul_pd(_Aelem, _Belem);
        _sum = _mm_add_pd(_sum, _mul);
    }
    double *result = (double *)malloc(sizeof(double) * 2);
    _mm_store_pd(result, _sum);
    sum = result[0] + result[1];
    return sum;
}

int main()
{
    const int N = 1024;

    double *A = (double *)malloc(sizeof(double) * N * N);
    double *B = (double *)malloc(sizeof(double) * N * N);
  
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

    double start_time = omp_get_wtime();
    double *C = (double *)malloc(sizeof(double) * N * N);    
    double *B_tr = (double *)malloc(sizeof(double) * N * N);
    
    transpose(B, B_tr, N);
    free(B);

    #pragma omp parallel for //firstprivate(A, B_tr, N)
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            C[i * N + j] = convolute(A + i * N, B_tr + j * N, N);
        }
    }

    free(B_tr);
    free(A);
    free(C);
    double runtime = omp_get_wtime() - start_time;

    printf("\nruntime was %4.4lf\n", runtime);


    
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

    // return 0;    
}