double convolute(double *A, double *B, int n){
    double sum = 0.0;
    for (int i = 0; i < n; i++){
        sum += A[i] * B[i];
    }
    return sum;
}