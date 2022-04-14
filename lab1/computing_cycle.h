int computing_cycle(int size, int rank, double* f, double* U0, double* Solution, double C, double tau, int K, int M);
int ind(int k, int m, int width);
int template(double *Solution, double *f, int i0, int i1, int i2, int i3, double C, double tau);
int template_left(double* Solution, double *f, int i1, int i2, int i3, double C, double tau);
int template_right(double* Solution, double *f, int i1, int i2, int i3, double C, double tau);
