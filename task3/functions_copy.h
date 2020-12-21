#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <cmath>
//#include <sys/time.h>



//double get_full_time();

double norm (const int n, double *array);
double f (int k, int n, int i, int j);
void fillMatrix (int n, double *array, int k, int total_p, int my_rank);
int readMatrix (int n, double *array, std::string filename, int total_p, int my_rank);

int fillFullMatrix (double *array, int n, int k, std::string filename);
void fillVector (double *vector, double *array, int n);

void print (double *array, int l, int n, int m);
void residual (const int n, double *array, double *vector, double *x);
void error (const int n, double *x);


#endif 
