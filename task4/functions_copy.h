#ifndef FUNCTIONS_H
#define FUNCTIONS_H



int fillMatrix (double *array, double *array_copy, int n, int k, std::string filename);
double euclidNorm (const int n, double *x);
double matrixNorm (const int n, double *array);
double f(int k, int n, int i, int j);

int reflect (int n, double *array, double *x, double *z);
int n_(int n, double *array, double alpha);  // n_(a)==S(A-alpha*E)
void recur (int n, double *array, double eps, double *spector, int *index, double a, double b);

void print (double *array, int l, int n, int m);
double error1 (const int n, double *array_copy, double *spector);
double error2 (const int n, double *array_copy, double *spector);

#endif 