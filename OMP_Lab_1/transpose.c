#include <stdlib.h>
int transpose(double *Matrix, double *Matrix_tr, int n){    
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            Matrix_tr[j * n + i] = Matrix[i * n + j];
        }        
    }
    return 0;
}