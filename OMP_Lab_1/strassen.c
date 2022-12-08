#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "print_matrix.h"
#include "transpose.h"
#include "convolute.h"
#include <omp.h>

int main()
{
    const int N = 1024;

    double *A = (double *)calloc(sizeof(double), N * N);
    double *B = (double *)calloc(sizeof(double), N * N);
  
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            A[i * N + j] = (i * N + j);
            if (i == j){
                B[i * N + j] = 10.0;
            } else {
                B[i * N + j] = 0.0;
            }
        }
    }

    double *C = (double *)calloc(sizeof(double), N * N);  

    double start_time = omp_get_wtime();
    strassen(A, B, C, N);
    double runtime = omp_get_wtime() - start_time;

    printf("\nruntime was %4.4lf\n", runtime);

    // print_matrix(A, N);
    // print_matrix(B, N);
    // print_matrix(C, N);

    free(A);
    free(C);
    free(B);    

    return 0;    
}

void matmul_2(double* A, double* B, double* C){
    double M1 = (A[0] + A[3]) * (B[0] + B[3]);
    double M2 = (A[2] + A[3]) * B[0];
    double M3 = A[0] * (B[1] - B[3]);
    double M4 = A[3] * (B[2] - B[0]);
    double M5 = (A[0] + A[1]) * B[3];
    double M6 = (A[2] - A[0]) * (B[0] + B[1]);
    double M7 = (A[1] - A[3]) * (B[2] + B[3]);
    C[0] = M1 + M4 - M5 + M7;
    C[1] = M3 + M5;
    C[2] = M2 + M4;
    C[3] = M1 - M2 + M3 + M6;
}

void add(double* A, double* B, double* C, int N){
    for (int i = 0; i < N * N; i++){
        C[i] = A[i] + B[i];
    }
}

void sub(double* A, double* B, double* C, int N){
    for (int i = 0; i < N * N; i++){
        C[i] = A[i] - B[i];
    }
}

void strassen(double *A, double *B, double *C, int N){

	int new_N = N/2;
	
    if(new_N == 1){
        matmul_2(A,B,C);
    }
    else {
        
        double* a11 = (double *)calloc(sizeof(double), new_N * new_N);
        double* a12 = (double *)calloc(sizeof(double), new_N * new_N);
        double* a21 = (double *)calloc(sizeof(double), new_N * new_N);
        double* a22 = (double *)calloc(sizeof(double), new_N * new_N);

        double* b11 = (double *)calloc(sizeof(double), new_N * new_N);
        double* b12 = (double *)calloc(sizeof(double), new_N * new_N);
        double* b21 = (double *)calloc(sizeof(double), new_N * new_N);
        double* b22 = (double *)calloc(sizeof(double), new_N * new_N);

        double* c11 = (double *)calloc(sizeof(double), new_N * new_N);
        double* c12 = (double *)calloc(sizeof(double), new_N * new_N);
        double* c21 = (double *)calloc(sizeof(double), new_N * new_N);
        double* c22 = (double *)calloc(sizeof(double), new_N * new_N);

        double* wAM1 = (double *)calloc(sizeof(double), new_N * new_N);
        double* wBM1 = (double *)calloc(sizeof(double), new_N * new_N);
        double* wAM2 = (double *)calloc(sizeof(double), new_N * new_N);
        double* wBM3 = (double *)calloc(sizeof(double), new_N * new_N);
        double* wBM4 = (double *)calloc(sizeof(double), new_N * new_N);
        double* wAM5 = (double *)calloc(sizeof(double), new_N * new_N);
        double* wAM6 = (double *)calloc(sizeof(double), new_N * new_N);
        double* wBM6 = (double *)calloc(sizeof(double), new_N * new_N);
        double* wAM7 = (double *)calloc(sizeof(double), new_N * new_N);
        double* wBM7 = (double *)calloc(sizeof(double), new_N * new_N);

        double* m1 = (double *)calloc(sizeof(double), new_N * new_N);
        double* m2 = (double *)calloc(sizeof(double), new_N * new_N);
        double* m3 = (double *)calloc(sizeof(double), new_N * new_N);
        double* m4 = (double *)calloc(sizeof(double), new_N * new_N);
        double* m5 = (double *)calloc(sizeof(double), new_N * new_N);
        double* m6 = (double *)calloc(sizeof(double), new_N * new_N);
        double* m7 = (double *)calloc(sizeof(double), new_N * new_N);


	    int i ,j;
	    for(i = 0;i < new_N; ++i)
	    {
	        for(j = 0; j < new_N;++j){
	            a11[i*new_N + j] = A[i*N + j];
	            a12[i*new_N + j] = A[i*N + j + new_N];
	            a21[i*new_N + j] = A[(i + new_N)*N + j];
	            a22[i*new_N + j] = A[(i + new_N)*N + j + new_N];

                b11[i*new_N + j] = B[i*N + j];
	            b12[i*new_N + j] = B[i*N + j + new_N];
	            b21[i*new_N + j] = B[(i + new_N)*N + j];
	            b22[i*new_N + j] = B[(i + new_N)*N + j + new_N];
	        }
	    }
        
        #pragma omp parallel 
        {
            #pragma omp single
            {
                #pragma omp task// M1 = (A11 + A22)*(B11 + B22)
                { 
                    add(a11, a22, wAM1, new_N);
                    add(b11, b22, wBM1, new_N);
                    strassen(wAM1,wBM1,m1, new_N);
                }
                #pragma omp task//M2 = (A21 + A22)*B11
                { 
                    add(a21, a22,wAM2, new_N);
                    strassen(wAM2,b11,m2, new_N);
                }
                #pragma omp task//M3 = A11*(B12 - B22)
                {
                    sub(b12, b22, wBM3, new_N);
                    strassen(a11, wBM3, m3, new_N);
                }
                #pragma omp task//M4 = A22*(B21 - B11)
                { 
                    sub(b21, b11, wBM4, new_N);
                    strassen(a22, wBM4, m4, new_N);
                }
                #pragma omp task//M5 = (A11 + A12)*B22
                {
                    add(a11, a12, wAM5, new_N);
                    strassen(wAM5, b22, m5, new_N);
                }
                #pragma omp task//M6 = (A21 - A11)*(B11 + B12)
                {
                    sub(a21, a11, wAM6, new_N);
                    add(b11, b12, wBM6, new_N);
                    strassen(wAM6, wBM6, m6, new_N);
                }
                #pragma omp task//M7 = (A12 - A22)*(B21 + B22)
                { 
                    sub(a12, a22, wAM7, new_N);
                    add(b21, b22, wBM7, new_N);
                    strassen(wAM7, wBM7, m7, new_N);
                }

                #pragma omp taskwait
                for(i = 0; i< new_N; i++){
                    for(j = 0; j<new_N; j++){
                        c11[i*new_N + j] = m1[i*new_N + j] + m4[i*new_N + j] - m5[i*new_N + j]+ m7[i*new_N + j];
                        c12[i*new_N + j] = m3[i*new_N + j] + m5[i*new_N + j];
                        c21[i*new_N + j] = m2[i*new_N + j] + m4[i*new_N + j];
                        c22[i*new_N + j] = m1[i*new_N + j] - m2[i*new_N + j] + m3[i*new_N + j] + m6[i*new_N + j];
                    }
                }

                for(i = 0; i < new_N; ++i){
                    for(j = 0; j < new_N; ++j){
                        C[i*N + j] = c11[i*new_N + j];
                        C[i*N + j + new_N] = c12[i*new_N + j];
                        C[(i + new_N)*N + j] = c21[i*new_N + j];
                        C[(i + new_N)*N + j + new_N] = c22[i*new_N + j];

                    }
                }
            }
        } 
        free (a11);
        free (a12);
        free (a21);
        free (a22);

        free (b11);
        free (b12);
        free (b21);
        free (b22);

        free (c11);
        free (c12);
        free (c21);
        free (c22);

        free(wAM1);
        free(wBM1);
        free(wAM2);
        free(wBM3);
        free(wBM4);
        free(wAM5);
        free(wAM6);
        free(wBM6);
        free(wAM7);
        free(wBM7);

        free(m1);
        free(m2);
        free(m3);
        free(m4);
        free(m5);
        free(m6);
        free(m7);
    }    
}