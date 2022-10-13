#include <stdlib.h>
double *reorder(double *Matrix, int n){
    double *result = (double *)malloc(sizeof(double) * n * n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            result[j * n + i] = Matrix[i * n + j];
        }        
    }
    return result;
}