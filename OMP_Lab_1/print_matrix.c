#include <stdio.h>
#include "print_matrix.h"

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