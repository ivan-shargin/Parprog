#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

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

float convolute(float *A, float *B, int n){
    float sum = 0.0;
    for (int i = 0; i < n; i++){
        sum += A[i] * B[i];
    }
    return sum;
}

int main()
{
    const int N = 1024;

    float *A = (float *)calloc(sizeof(float), N * N);
    float *B = (float *)calloc(sizeof(float), N * N);
  
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

    float start_time = omp_get_wtime();
    float *C = (float *)calloc(sizeof(float), N * N);    
    float *B_tr = (float *)calloc(sizeof(float), N * N);
    
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
    float runtime = omp_get_wtime() - start_time;

    printf("\nruntime was %4.4lf\n", runtime);


    
    // float A_test[4] = {10.0, 0.0, 0.0, 10.0};
    // float B_test[4] = {1.0, 2.0, 3.0, 4.0};
    // int n = 2;
    // print_matrix(A_test, n);
    // print_matrix(B_test, n);

    // float C_test[4];
    // float *B_test_tr = (float *)malloc(sizeof(double) * n * n);
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