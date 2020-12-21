#include <iostream>
#include <string>
#include <cmath>
#include "functions_copy.h"



double get_full_time()
{
	struct timeval r;
	gettimeofday(&r,0);
	return (double)r.tv_sec+(double)r.tv_usec/1000000.;
}

double norm (const int n, double *array) // standard norm in R^(n*n)
{
    double value = 0;
    for (int i=0; i<n*n; i++)
    {
        value += array[i]*array[i];
    }
    return sqrt(value);
}


double f(int k, int n, int i, int j)  // i,j = 1, ..., n (NOT from 0!)
{
    switch (k)
    {
        case 1:
            return n - ((i<j)?j:i) + 1;
        case 2:
            return (i<j)?j:i;
        case 3:
            return (i<j)?(j-i):(i-j);
        case 4:
            return 1./(i+j-1);
    }
    return 0;
}


int fillMatrix (double *array, int n, int k, std::string filename)
{
    if (k == 0)
    {
        FILE *fp = fopen (filename.c_str(), "r");
        if (!fp)
            {
                std::cout << "File cannot be opened\n";
                return -1;
            }
        int i;
        for (i=0; i < n*n; i++)
        {
            if (!fscanf (fp, "%lf", array + i))
            {
                std::cout << "File cannot be read / data is invalid\n";
                return -1;
            }
        }
        fclose(fp);
    }
    else
    {
        for (int i=0; i < n*n; i++)
        {
            array[i] = f(k,n,(i/n)+1,(i%n)+1);
        }
    }
    return 0;
}


void fillVector (double *vector, double *b, double *array, int n)
{
    double temp;
    for (int i=0; i<n; i++)
    {
        temp = 0;
        for (int k=0; k <= (n-1)/2; k++)
        {
            temp += array[i*n + 2*k];     // a(i+1)(2k+1)  --->  array[i*n + 2k]
        }
        vector[i] = temp;
        b[i] = temp;
    }
}


void print (double *array, int l, int n, int m)
{
    for (int i=0; i < ((l<m)?l:m); i++)
    {
        for (int j=0; j < ((n<m)?n:m); j++)
        {
            printf(" %10.3e", array[i*n + j]); 
        }
        printf("\n");
    }
    printf("\n");
}


void residual (const int n, double *A, double *b, double *x)
{
    double value;

    double *Ax = new (std::nothrow) double[n];
    if (!Ax)
    {
        std::cout << "couldn't allocate memory\n";
        return;
    }
    for (int i=0; i<n; i++)
    {
        value = 0;
        for (int j=0; j<n; j++)
        {
            value += A[i*n + j]*x[j];
        }
        Ax[i] = value;
    }
    double norm1 = 0, norm2 = 0; // norm1 = ||Ax-b||^2     norm2 = ||b||^2
    for (int i=0; i<n; i++)
    {
        norm1 += (Ax[i]-b[i])*(Ax[i]-b[i]);
        norm2 += b[i]*b[i];
    }

    std::cout << " Residual (||Ax-b|| / ||b||) of the algorithm: " << sqrt(norm1)/sqrt(norm2) << "\n";

    delete[] Ax;

}


void error (const int n, double *x)
{
    double value = 0;

    for (int i=0; i<n; i++)
    {
        value += ((i+1)%2 - x[i])*((i+1)%2 - x[i]);
    }
    std::cout << " Error norm is: " << sqrt(value) << "\n";
}


