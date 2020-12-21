#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <pthread.h>
#include <cmath>
#include <sys/time.h>


class Arg 
{
public:

	int n;  				// matrix size
	int p;  				// total number of threads
	int my_rank;  			// number of current thread
	int *err;  				// errors 
	double *A;  			// pointer to matrix
	double *b;  			// pointer to vector 
	double *x;  			// pointer to solution 
	double *cosphi;			// common value
	double *sinphi;			// common value
	double time_thr;  		// time_thr - execution time for current thread 
	pthread_barrier_t *barrier;
	//pthread_mutex_t *mutex;
	int *flag;
	
	Arg ()
	{
		n = 0;
		p = 0;
		my_rank = 0;
		err = 0;
		A = nullptr;
		b = nullptr;
		x = nullptr;
		cosphi = nullptr;
		sinphi = nullptr;
		time_thr = 0;
		flag = nullptr;
	}
};


//double get_time();
double get_full_time();


double f (int k, int n, int i, int j);
// array - initial matrix, A - its copy
int fillMatrix (double *array, int n, int k, std::string filename);
// vector - initial b, b - its copy
void fillVector (double *vector, double *b, double *array, int n);
void print (double *array, int l, int n, int m);
double norm (const int n, double *array);
void residual (const int n, double *A, double *vector, double *x);
void error (const int n, double *x);


#endif 
