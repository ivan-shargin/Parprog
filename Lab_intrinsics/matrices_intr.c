#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <emmintrin.h>
float result[4];

int print_matrix(float *Matrix, int n){
    printf("\n");
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            printf("%8.1lf", Matrix[i * n + j]);
        }
        printf("\n");
    }
    return 0;
}

int transpose(float *Matrix, float *Matrix_tr, int n){    
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            Matrix_tr[j * n + i] = Matrix[i * n + j];
        }        
    }
    return 0;
}

float convolute(float *A, float *B, int n, float *result){
    __m128 _sum, _Aelem, _Belem, _mul;
    _sum = _mm_set_ps(0.0, 0.0, 0.0, 0.0);
    for (int i = 0; i < n; i+=4){
        _Aelem = _mm_set_ps(A[i], A[i+1], A[i+2], A[i+3]);
        _Belem = _mm_set_ps(B[i], B[i+1], B[i+2], B[i+3]);
        _mul = _mm_mul_ps(_Aelem, _Belem);
        _sum = _mm_add_ps(_sum, _mul);
    }
    
    _mm_store_ps(result, _sum);
    float sum = result[0] + result[1] + result[2] + result[3];
    return sum;
}

int main()
{
    const int N = 1024;

    float *A = (float *)aligned_alloc(sizeof(float),sizeof(float) * N * N);
    float *B = (float *)aligned_alloc(sizeof(float),sizeof(float) * N * N);
  
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

    // print_matrix(A, N);
    // print_matrix(B, N);

    float start_time = omp_get_wtime();
    float *C = (float *)aligned_alloc(sizeof(float),sizeof(float) * N * N);    
    float *B_tr = (float *)aligned_alloc(sizeof(float),sizeof(float) * N * N);
    float *result = (float*)malloc(sizeof(float)*4);
    
    transpose(B, B_tr, N);
    free(B);

    #pragma omp parallel for //firstprivate(A, B_tr, N)
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            C[i * N + j] = convolute(A + i * N, B_tr + j * N, N, result);
        }
    }

   
    // print_matrix(C, N);

    free(B_tr);
    free(A);
    free(C);
    free(result);
    float runtime = omp_get_wtime() - start_time;

    printf("\nruntime was %4.4lf\n", runtime);

    


    
    // float A_test[16] = {10.0, 0.0,0. 0.0, 10.0};
    // float B_test[16] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0};
    // int n = 4;
    // print_matrix(A_test, n);
    // print_matrix(B_test, n);

    // float C_test[4];
    // float *B_test_tr = (float *)aligned_alloc(sizeof(float),sizeof(float) * N * N);
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