#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "print_matrix.h"
#include "transpose.h"
#include "convolute.h"

int main()
{
    double A[4] = {10.0, 0.0, 0.0, 10.0};
    double B[4] = {1.0, 2.0, 3.0, 4.0};
    int n = 2;
    print_matrix(A, n);
    // print_matrix(transpose(B, n), n);
    double C[4];
    double *B_tr = (double *)malloc(sizeof(double) * n * n);
    transpose(B, B_tr, n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            C[i * n + j] = convolute(A + i * n, B_tr + j * n, n);
        }
    }

    print_matrix(C, n);
    free(B_tr);

    return 0;    
}