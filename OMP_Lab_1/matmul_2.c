void matmul_2(double* A, double* B, double* C){
    double M1 = (A[0] + A[3]) * (B[0] + B[3]);
    double M2 = (A[2] + A[3] * B[0]);
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